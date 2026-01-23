!> A tracer package for passive dye tracers with constant surface fluxes in specified regions.
module dye_fluxes

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_coms,               only : EFP_type
use MOM_coupler_types,      only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator,      only : diag_ctrl, post_data, register_diag_field
use MOM_error_handler,      only : MOM_error, FATAL, WARNING
use MOM_file_parser,        only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,       only : forcing
use MOM_grid,               only : ocean_grid_type
use MOM_hor_index,          only : hor_index_type
use MOM_interface_heights,  only : thickness_to_dz
use MOM_io,                 only : vardesc, var_desc, query_vardesc
use MOM_open_boundary,      only : ocean_OBC_type
use MOM_restart,            only : query_initialized, MOM_restart_CS
use MOM_spatial_means,      only : global_mass_int_EFP
use MOM_sponge,             only : set_up_sponge_field, sponge_CS
use MOM_time_manager,       only : time_type
use MOM_tracer_registry,    only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic,    only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init,      only : tracer_Z_init
use MOM_unit_scaling,       only : unit_scale_type
use MOM_variables,          only : surface, thermo_var_ptrs
use MOM_verticalGrid,       only : verticalGrid_type
use MOM_tracer_advect_schemes, only : set_tracer_advect_scheme, TracerAdvectionSchemeDoc

implicit none ; private

#include <MOM_memory.h>

public register_dye_flux_tracer, initialize_dye_flux_tracer
public dye_flux_tracer_column_physics, dye_flux_tracer_surface_state
public dye_flux_stock, dye_fluxes_end

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the dye flux tracer package
type, public :: dye_flux_tracer_CS ; private
  integer :: ntr    !< The number of tracers that are actually used.
  logical :: coupled_tracers = .false.  !< These tracers are not offered to the coupler.
  
  real, allocatable, dimension(:) :: flux_source_minlon !< Minimum longitude of region where flux is applied
                                                        !! [degrees_E] or [km] or [m]
  real, allocatable, dimension(:) :: flux_source_maxlon !< Maximum longitude of region where flux is applied
                                                        !! [degrees_E] or [km] or [m]
  real, allocatable, dimension(:) :: flux_source_minlat !< Minimum latitude of region where flux is applied
                                                        !! [degrees_N] or [km] or [m]
  real, allocatable, dimension(:) :: flux_source_maxlat !< Maximum latitude of region where flux is applied
                                                        !! [degrees_N] or [km] or [m]
  real, allocatable, dimension(:) :: surface_flux_const !< Constant surface flux value for each tracer
                                                        !! [conc Z T-1 ~> conc m s-1]
  
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the tracer registry
  real, pointer :: tr(:,:,:,:) => NULL() !< The array of tracers used in this subroutine [CU ~> conc]

  integer, allocatable, dimension(:) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux if it is used and the
                                               !! surface tracer concentrations are to be provided to the coupler.

  integer, allocatable, dimension(:) :: id_surface_flux !< Diagnostic IDs for surface tracer fluxes
  integer, allocatable, dimension(:) :: id_tr_dia_diff  !< Diagnostic IDs for vertical tracer fluxes (positive up)

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc), allocatable :: tr_desc(:) !< Descriptions and metadata for the tracers
  logical :: tracers_may_reinit = .true. !< If true the tracers may be initialized if not found in a restart file
end type dye_flux_tracer_CS

contains

!> This subroutine is used to register tracer fields and subroutines
!! to be used with MOM.
function register_dye_flux_tracer(HI, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure.
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),      intent(in) :: US   !< A dimensional unit scaling type
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(dye_flux_tracer_CS),   pointer    :: CS   !< A pointer that is set to point to the control
                                                 !! structure for this module
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                 !! structure for the tracer advection and diffusion module.
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control structure

  ! Local variables
  character(len=40)  :: mdl = "dye_fluxes" ! This module's name.
  character(len=48)  :: var_name ! The variable's name.
  character(len=48)  :: desc_name ! The variable's descriptor.
  character(len=48)  :: param_name ! The param's name suffix.
  character(len=48)  :: flux_units ! The units for tracer fluxes
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  real, pointer :: tr_ptr(:,:,:) => NULL() ! A pointer to one of the tracers [CU ~> conc]
  logical :: register_dye_flux_tracer
  integer :: isd, ied, jsd, jed, nz, m
  integer :: advect_scheme   ! Advection scheme value for this tracer
  character(len=256) :: mesg ! Advection scheme name for this tracer

  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(FATAL, "register_dye_flux_tracer called with an "// &
                          "associated control structure.")
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NUM_DYE_FLUX_TRACERS", CS%ntr, &
                 "The number of dye tracers with surface fluxes in this run. "//&
                 "Each tracer should have a separate surface flux region.", default=0)
  
  allocate(CS%flux_source_minlon(CS%ntr), &
           CS%flux_source_maxlon(CS%ntr), &
           CS%flux_source_minlat(CS%ntr), &
           CS%flux_source_maxlat(CS%ntr), &
           CS%surface_flux_const(CS%ntr))
  allocate(CS%ind_tr(CS%ntr))
  allocate(CS%tr_desc(CS%ntr))
  allocate(CS%id_tr_dia_diff(CS%ntr))
  allocate(CS%id_surface_flux(CS%ntr))

  ! Read geographic bounds and flux values for each tracer
  call get_param(param_file, mdl, "DYE_FLUX_SOURCE_MINLON", CS%flux_source_minlon, &
                 "This is the minimum longitude of the region where the dye flux is applied.", &
                 units=G%x_ax_unit_short, fail_if_missing=.true.)
  if (minval(CS%flux_source_minlon(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_flux_tracer: Not enough values provided for DYE_FLUX_SOURCE_MINLON")

  call get_param(param_file, mdl, "DYE_FLUX_SOURCE_MAXLON", CS%flux_source_maxlon, &
                 "This is the maximum longitude of the region where the dye flux is applied.", &
                 units=G%x_ax_unit_short, fail_if_missing=.true.)
  if (minval(CS%flux_source_maxlon(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_flux_tracer: Not enough values provided for DYE_FLUX_SOURCE_MAXLON")

  call get_param(param_file, mdl, "DYE_FLUX_SOURCE_MINLAT", CS%flux_source_minlat, &
                 "This is the minimum latitude of the region where the dye flux is applied.", &
                 units=G%y_ax_unit_short, fail_if_missing=.true.)
  if (minval(CS%flux_source_minlat(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_flux_tracer: Not enough values provided for DYE_FLUX_SOURCE_MINLAT")

  call get_param(param_file, mdl, "DYE_FLUX_SOURCE_MAXLAT", CS%flux_source_maxlat, &
                 "This is the maximum latitude of the region where the dye flux is applied.", &
                 units=G%y_ax_unit_short, fail_if_missing=.true.)
  if (minval(CS%flux_source_maxlat(:)) < -1.e29) &
    call MOM_error(FATAL, "register_dye_flux_tracer: Not enough values provided for DYE_FLUX_SOURCE_MAXLAT")

  call get_param(param_file, mdl, "DYE_FLUX_SURFACE_VALUE", CS%surface_flux_const, &
                 "This is the constant surface flux value for each dye tracer. "//&
                 "Positive values represent flux into the ocean.", &
                 units="conc m s-1", scale=US%m_to_Z*US%s_to_T, fail_if_missing=.true.)
  if (minval(CS%surface_flux_const(:)) < -1.e29*US%m_to_Z*US%s_to_T) &
    call MOM_error(FATAL, "register_dye_flux_tracer: Not enough values provided for DYE_FLUX_SURFACE_VALUE")

  allocate(CS%tr(isd:ied,jsd:jed,nz,CS%ntr), source=0.0)

  do m = 1, CS%ntr
    write(param_name(:),'(A,I3.3,A)') "DYE_FLUX",m,"_TRACER_ADVECTION_SCHEME"
    call get_param(param_file, mdl, trim(param_name), mesg, &
          desc="The horizontal transport scheme for dye flux tracer:\n"//&
          trim(TracerAdvectionSchemeDoc)//&
          "\n Set to blank (the default) to use TRACER_ADVECTION_SCHEME.", default="")
    call set_tracer_advect_scheme(mesg, advect_scheme)

    ! Register each tracer
    write(var_name,'(A,I3.3)') "dye_flux_",m
    write(desc_name,'(A,I3.3)') "Concentration of Dye Flux Tracer ",m
    CS%tr_desc(m) = var_desc(var_name, units="kg kg-1", longname=desc_name, caller=mdl)

    ! This needs to be changed if the units of tracer are changed above.
    if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
    else ; flux_units = "kg s-1" ; endif

    tr_ptr => CS%tr(:,:,:,m)
    call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, &
                         tr_desc=CS%tr_desc(m), registry_diags=.true., &
                         restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit, &
                         advect_scheme=advect_scheme)

    !   Set coupled_tracers to be true (hard-coded above) to provide the surface
    ! values to the coupler (if any).  This is meta-code and its arguments will
    ! currently (deliberately) give fatal errors if it is used.
    if (CS%coupled_tracers) &
      CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(var_name)//'_flux', &
          flux_type=' ', implementation=' ', caller="register_dye_flux_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_dye_flux_tracer = .true.
end function register_dye_flux_tracer

!> This subroutine initializes the CS%ntr tracer fields in tr(:,:,:,:)
!! and it sets up the tracer output.
subroutine initialize_dye_flux_tracer(restart, day, G, GV, US, h, diag, OBC, CS, sponge_CSp, tv)
  logical,                            intent(in) :: restart !< .true. if the fields have already been
                                                            !! read from a restart file.
  type(time_type), target,            intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  type(unit_scale_type),              intent(in) :: US   !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl), target,            intent(in) :: diag !< Structure used to regulate diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(dye_flux_tracer_CS),           pointer    :: CS   !< The control structure returned by a previous
                                                         !! call to register_dye_flux_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp !< A pointer to the control structure
                                                         !! for the sponges, if they are in use.
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure pointing to various thermodynamic variables

  ! Local variables
  character(len=64)  :: var_name, longname
  integer :: i, j, k, m

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  CS%diag => diag

  ! Register diagnostics for surface flux and vertical diffusive flux
  do m = 1, CS%ntr
    write(var_name,'(A,I3.3,A)') "dye_flux",m,"_sflux"
    write(longname,'(A,I3.3,A)') "Surface flux of dye flux tracer ",m," (positive into ocean)"
    CS%id_surface_flux(m) = register_diag_field('ocean_model', trim(var_name), &
        diag%axesT1, day, trim(longname), 'conc m s-1', conversion=US%Z_to_m*US%s_to_T)

    write(var_name,'(A,I3.3,A)') "dye_flux",m,"_dia_diff"
    write(longname,'(A,I3.3,A)') "Vertical diffusive flux of dye flux tracer ",m," (positive up)"
    CS%id_tr_dia_diff(m) = register_diag_field('ocean_model', trim(var_name), &
        diag%axesTi, day, trim(longname), 'conc H s-1', conversion=GV%H_to_MKS*US%s_to_T)
  enddo

  ! Initialize tracers to zero (can be modified if restart file exists)
  if (.not. restart) then
    do m = 1, CS%ntr
      do k = 1, GV%ke ; do j = G%jsc, G%jec ; do i = G%isc, G%iec
        CS%tr(i,j,k,m) = 0.0
      enddo ; enddo ; enddo
    enddo
  endif

end subroutine initialize_dye_flux_tracer

!> This subroutine applies diapycnal diffusion, surface fluxes, and any other column
!! tracer physics or chemistry to the tracers from this file.
subroutine dye_flux_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, tv, CS, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  type(dye_flux_tracer_CS), pointer   :: CS   !< The control structure returned by a previous
                                              !! call to register_dye_flux_tracer.
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work ! Used so that h can be modified [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)+1) :: vert_flux ! Vertical tracer flux positive upward
                                              !! [conc H T-1 ~> conc m s-1]
  real, dimension(SZI_(G),SZJ_(G)) :: surface_flux ! Surface tracer flux [conc Z T-1 ~> conc m s-1]
  real    :: Idt      ! Inverse of timestep [T-1 ~> s-1]
  real    :: flux_mag ! Magnitude of surface flux [conc Z T-1 ~> conc m s-1]
  integer :: i, j, k, is, ie, js, je, nz, m

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  Idt = 1.0 / dt

  ! Apply surface fluxes and vertical diffusion for each tracer
  do m = 1, CS%ntr
    ! Initialize surface flux field
    surface_flux(:,:) = 0.0
    
    ! Apply constant surface flux in the specified region
    do j = js, je ; do i = is, ie
      if (CS%flux_source_minlon(m) < G%geoLonT(i,j) .and. &
          CS%flux_source_maxlon(m) >= G%geoLonT(i,j) .and. &
          CS%flux_source_minlat(m) < G%geoLatT(i,j) .and. &
          CS%flux_source_maxlat(m) >= G%geoLatT(i,j) .and. &
          G%mask2dT(i,j) > 0.0 ) then
        surface_flux(i,j) = CS%surface_flux_const(m)
      endif
    enddo ; enddo

    ! Post diagnostic of surface flux
    if (CS%id_surface_flux(m) > 0) &
      call post_data(CS%id_surface_flux(m), surface_flux, CS%diag)

    ! Apply surface flux to top layer
    ! Surface flux modifies tracer concentration: dC/dt = Flux / h
    do j = js, je ; do i = is, ie
      if (G%mask2dT(i,j) > 0.0 .and. h_old(i,j,1) > 0.0) then
        CS%tr(i,j,1,m) = CS%tr(i,j,1,m) + (dt * surface_flux(i,j) / h_old(i,j,1))
      endif
    enddo ; enddo

    ! Apply vertical diffusion
    if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
      do k = 1, nz ; do j = js, je ; do i = is, ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
                                          evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    else
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    endif

    ! Calculate net vertical flux from entrainment for diagnostics
    ! Net flux = upward component - downward component
    ! Upward (from below): eb(k) * tr(k+1), Downward (from above): ea(k+1) * tr(k)
    do K = 2, nz ; do j = js, je ; do i = is, ie
      vert_flux(i,j,K) = 0.0
      if (k < nz) vert_flux(i,j,K) = vert_flux(i,j,k) + (eb(i,j,k) * CS%tr(i,j,k+1,m)) * Idt
      if (k < nz) vert_flux(i,j,K) = vert_flux(i,j,k) - (ea(i,j,k+1) * CS%tr(i,j,k,m)) * Idt
    enddo ; enddo ; enddo
    do j = js, je ; do i = is, ie ; vert_flux(i,j,1) = 0.0 ; vert_flux(i,j,nz+1) = 0.0 ; enddo ; enddo

    ! Post diagnostic of vertical flux
    if (CS%id_tr_dia_diff(m) > 0) &
      call post_data(CS%id_tr_dia_diff(m), vert_flux, CS%diag)
  enddo

end subroutine dye_flux_tracer_column_physics

!> This function calculates the mass-weighted integral of all tracer stocks,
!! returning the number of stocks it has calculated.  If the stock_index
!! is present, only the stock corresponding to that coded index is returned.
function dye_flux_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),              intent(in)    :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in)    :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),       intent(out)   :: stocks !< The mass-weighted integrated amount of each
                                                            !! tracer, in kg times concentration units [kg conc]
  type(dye_flux_tracer_CS),           pointer       :: CS   !< The control structure returned by a
                                                            !! previous call to register_dye_flux_tracer.
  character(len=*), dimension(:),     intent(out)   :: names  !< The names of the stocks calculated.
  character(len=*), dimension(:),     intent(out)   :: units  !< The units of the stocks calculated.
  integer, optional,                  intent(in)    :: stock_index !< The coded index of a specific stock
                                                            !! being sought.
  integer :: dye_flux_stock !< The number of stocks calculated here.

  ! Local variables
  integer :: m

  dye_flux_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    ! Check whether this stock is available from this routine.

    ! No stocks from this routine are being checked yet.  Return 0.
    return
  endif ; endif

  do m = 1, CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), caller="dye_flux_stock")
    units(m) = trim(units(m))//" kg"
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  enddo
  dye_flux_stock = CS%ntr

end function dye_flux_stock

!> This subroutine extracts the surface fields from this tracer package that
!! are to be shared with the atmosphere in coupled configurations.
!! This particular tracer package does not report anything back to the coupler.
subroutine dye_flux_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< A structure containing fields that
                                               !! describe the surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h  !< Layer thickness [H ~> m or kg m-2].
  type(dye_flux_tracer_CS), pointer      :: CS !< The control structure returned by a previous
                                               !! call to register_dye_flux_tracer.
  ! Local variables
  integer :: m, is, ie, js, je, isd, ied, jsd, jed

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m = 1, CS%ntr
      !   This call loads the surface values into the appropriate array in the
      ! coupler-type structure.
      call set_coupler_type_data(CS%tr(:,:,1,m), CS%ind_tr(m), sfc_state%tr_fields, &
                   idim=(/isd, is, ie, ied/), jdim=(/jsd, js, je, jed/), turns=G%HI%turns)
    enddo
  endif

end subroutine dye_flux_tracer_surface_state

!> Clean up any allocated memory after the run.
subroutine dye_fluxes_end(CS)
  type(dye_flux_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                          !! call to register_dye_flux_tracer.

  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine dye_fluxes_end

!> \namespace dye_fluxes
!!
!!    This module contains a tracer package for passive dye tracers with
!!    constant surface fluxes applied over specified geographic regions.
!!    The surface flux is applied as a boundary condition in the column
!!    physics routine, modifying the tracer concentration in the top layer.
!!
!!    Key features:
!!    - Apply constant surface fluxes in user-defined lat/lon regions
!!    - Multiple independent dye tracers supported
!!    - Diagnostic output for surface and vertical fluxes
!!
!!    A single subroutine is called from within each file to register
!!    each of the tracers for reinitialization and advection and to
!!    register the subroutine that initializes the tracers and sets up
!!    their output and the subroutine that does any tracer physics or
!!    chemistry along with diapycnal mixing.

end module dye_fluxes
