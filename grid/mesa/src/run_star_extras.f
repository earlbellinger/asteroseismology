! ***********************************************************************
!
!   Copyright (C) 2010  Bill Paxton
!
!   this file is part of mesa.
!
!   mesa is free software; you can redistribute it and/or modify
!   it under the terms of the gnu general library public license as published
!   by the free software foundation; either version 2 of the license, or
!   (at your option) any later version.
!
!   mesa is distributed in the hope that it will be useful, 
!   but without any warranty; without even the implied warranty of
!   merchantability or fitness for a particular purpose.  see the
!   gnu library general public license for more details.
!
!   you should have received a copy of the gnu library general public license
!   along with this software; if not, write to the free software
!   foundation, inc., 59 temple place, suite 330, boston, ma 02111-1307 usa
!
! ***********************************************************************
 
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib
      
      ! real(dp) :: log_L_prev
      ! real(dp) :: logT_prev
      ! real(dp) :: d_log_L_dT
      
      implicit none
      
      ! these routines are called by the standard run_star check_model
      contains
      
       subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

         ! Uncomment these lines if you wish to use the functions in this file,
         ! otherwise we use a null_ version which does nothing.
         s% extras_startup => extras_startup
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         ! Once you have set the function pointers you want,
         ! then uncomment this (or set it in your star_job inlist)
         ! to disable the printed warning message,
          s% job% warn_run_star_extras =.false.

      end subroutine extras_controls

      ! None of the following functions are called unless you set their
      ! function point in extras_control.


      integer function extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_startup = 0
         if (.not. restart) then
            call alloc_extra_info(s)
         else ! it is a restart
            call unpack_extra_info(s)
         end if
      end function extras_startup


      ! returns either keep_going, retry, backup, or terminate.
      integer function extras_check_model(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         
         real(dp) :: Lnuc_div_L
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going
         if (.false. .and. s% star_mass_h1 < 0.35d0) then
            ! stop when star hydrogen mass drops to specified level
            extras_check_model = terminate
            write(*, *) 'have reached desired hydrogen mass'
            return
         end if
         
         ! for finding ZAMS 
         Lnuc_div_L = s% L_nuc_burn_total / s% L_phot
         if (s% x_logical_ctrl(2) .and. Lnuc_div_L > s% Lnuc_div_L_zams_limit) then
             if (Lnuc_div_L < s% Lnuc_div_L_zams_limit + 0.0001) then
                 extras_check_model = terminate
                 write(*, *) 'reached ZAMS'
             else 
                 extras_check_model = backup
                 write(*, *) 'zeroing in on ZAMS'
             end if 
         end if

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizeable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
      end function extras_check_model


      integer function how_many_extra_history_columns(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 3
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, id_extra, n, names, vals, ierr)
         integer, intent(in) :: id, id_extra, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         
         real(dp) :: d_log_L_dT, d_log_L_dT_cutoff, d_log_L_dT_tol
         real(dp) :: log_L, log_L_old, log_Teff, log_Teff_old
         
         real(dp) :: h_exh_core_mass, h_exh_core_radius
         integer :: k
         
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !note: do NOT add the extras names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.
         
         log_L = safe_log10_cr(s% L_phot)
         log_L_old = safe_log10_cr(s% L_phot_old)
         log_Teff = safe_log10_cr(s% Teff)
         log_Teff_old = safe_log10_cr(s% Teff_old)
         
         d_log_L_dT = (log_L - log_L_old) / (log_Teff - log_Teff_old)
         
         names(1) = "d_log_L_dT"
         vals(1) = d_log_L_dT
         
         
         h_exh_core_mass = 0
         h_exh_core_radius = 0
         if (s% center_h1 < s% xa_central_lower_limit(1)) then
             do k = 1, s% nz
                ! write(*,*) "xa1k", s% xa(1,k)
                if (s% xa(1,k) <= s% xa_central_lower_limit(1)) then
                    ! write(*,*) "h_exh_core_mass = ", s% m(k) / ( s% star_mass * msol )
                    ! write(*,*) "h_exh_core_radius = ", s% r(k) / ( s% photosphere_r * rsol )
                    h_exh_core_mass = s% m(k) / ( s% star_mass * msol )
                    h_exh_core_radius = s% r(k) / ( s% photosphere_r * rsol )
                    exit
                end if
             end do
         end if
         
         names(2) = "h_exh_core_mass"
         vals(2) = h_exh_core_mass 
         names(3) = "h_exh_core_radius"
         vals(3) = h_exh_core_radius
         
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id, id_extra)
         use star_def, only: star_info
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 0
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, id_extra, n, nz, names, vals, ierr)
         use star_def, only: star_info, maxlen_profile_column_name
         use const_def, only: dp
         integer, intent(in) :: id, id_extra, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         !if (n /= 1) stop 'data_for_extra_profile_columns'
         !names(1) = 'beta'
         !do k = 1, nz
         !   vals(k,1) = s% Pgas(k)/s% P(k)
         !end do

      end subroutine data_for_extra_profile_columns


      ! returns either keep_going or terminate.
      ! note: cannot request retry or backup; extras_check_model can do that.
      integer function extras_finish_step(id, id_extra)
         integer, intent(in) :: id, id_extra
         integer :: ierr
         type (star_info), pointer :: s
         
         real(dp) :: d_log_L_dT, log_L, log_L_old, log_Teff, log_Teff_old
         
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_finish_step = keep_going
         call store_extra_info(s)

         ! to save a profile,
            ! s% need_to_save_profiles_now = .true.
         ! to update the star log,
            ! s% need_to_update_history_now = .true.
         
         ! write(*, *) s% diffusion_class_factor(1)
         if (s% x_ctrl(2) > 0 .and. s% do_element_diffusion) then
             if (s% diffusion_class_factor(1) > 0) then
                 ! write(*, *) 'decreasing'
                 s% diffusion_class_factor(:) = &
                     s% diffusion_class_factor(:) - s% x_ctrl(2)
             end if 
             if (s% diffusion_class_factor(1) <= 0) then 
                 s% do_element_diffusion = .false.
                 s% diffusion_class_factor(:) = 0.
                 ! write(*, *) 'turning off'
             end if
         end if
         
         log_L = safe_log10_cr(s% L_phot)
         log_L_old = safe_log10_cr(s% L_phot_old)
         log_Teff = safe_log10_cr(s% Teff)
         log_Teff_old = safe_log10_cr(s% Teff_old)
         
         d_log_L_dT = (log_L - log_L_old) / (log_Teff - log_Teff_old)
         if (s% x_ctrl(1) > -1d99 .and. d_log_L_dT > s% x_ctrl(1) &
             .and. log_Teff > 3.675 .and. log_Teff < 3.74) then
             extras_finish_step = terminate
             write(*,*) "stopping: d_log_L_dT > x_ctrl(1)"
         end if
         
         if (s% x_ctrl(3) > -1d99 .and. s% delta_Pg > 0 &
             .and. s% delta_Pg < s% x_ctrl(3)) then
             extras_finish_step = terminate
             write(*,*) "stopping: delta_Pg < x_ctrl(3)"
         end if
         
         if (s% x_logical_ctrl(1) .and. log_Teff > log_Teff_old) then
             extras_finish_step = terminate
             write(*,*) "stopping: log_Teff increased and x_logical_ctrl(1)"
         end if
         
         ! set x_integer_ctrl(1) = 4 to stop when the He flash begins
         ! set x_integer_ctrl(1) = 6 to stop when He is burning
         if (s% x_integer_ctrl(1) > -1 .and. &
                 s% phase_of_evolution >= s% x_integer_ctrl(1)) then
             extras_finish_step = terminate
             write(*,*) "stopping: phase_of_evolution >= x_integer_ctrl(1)"
         end if
         ! phase_starting = 0
         ! phase_early_main_seq = 1
         ! phase_mid_main_seq = 2
         ! phase_wait_for_he = 3
         ! phase_he_ignition_over = 4
         ! phase_he_igniting = 5
         ! phase_helium_burning = 6
         ! phase_carbon_burning = 7
         
         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
      end function extras_finish_step


      subroutine extras_after_evolve(id, id_extra, ierr)
         integer, intent(in) :: id, id_extra
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve


      ! routines for saving and restoring extra data so can do restarts

         ! put these defs at the top and delete from the following routines
         !integer, parameter :: extra_info_alloc = 1
         !integer, parameter :: extra_info_get = 2
         !integer, parameter :: extra_info_put = 3


      subroutine alloc_extra_info(s)
         integer, parameter :: extra_info_alloc = 1
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_alloc)
      end subroutine alloc_extra_info


      subroutine unpack_extra_info(s)
         integer, parameter :: extra_info_get = 2
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_get)
      end subroutine unpack_extra_info


      subroutine store_extra_info(s)
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         call move_extra_info(s,extra_info_put)
      end subroutine store_extra_info


      subroutine move_extra_info(s,op)
         integer, parameter :: extra_info_alloc = 1
         integer, parameter :: extra_info_get = 2
         integer, parameter :: extra_info_put = 3
         type (star_info), pointer :: s
         integer, intent(in) :: op

         integer :: i, j, num_ints, num_dbls, ierr

         i = 0
         ! call move_int or move_flg
         num_ints = i

         i = 0
         ! call move_dbl

         num_dbls = i

         if (op /= extra_info_alloc) return
         if (num_ints == 0 .and. num_dbls == 0) return

         ierr = 0
         call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
         if (ierr /= 0) then
            write(*,*) 'failed in star_alloc_extras'
            write(*,*) 'alloc_extras num_ints', num_ints
            write(*,*) 'alloc_extras num_dbls', num_dbls
            stop 1
         end if

         contains

         subroutine move_dbl(dbl)
            real(dp) :: dbl
            i = i+1
            select case (op)
            case (extra_info_get)
               dbl = s% extra_work(i)
            case (extra_info_put)
               s% extra_work(i) = dbl
            end select
         end subroutine move_dbl

         subroutine move_int(int)
            integer :: int
            i = i+1
            select case (op)
            case (extra_info_get)
               int = s% extra_iwork(i)
            case (extra_info_put)
               s% extra_iwork(i) = int
            end select
         end subroutine move_int

         subroutine move_flg(flg)
            logical :: flg
            i = i+1
            select case (op)
            case (extra_info_get)
               flg = (s% extra_iwork(i) /= 0)
            case (extra_info_put)
               if (flg) then
                  s% extra_iwork(i) = 1
               else
                  s% extra_iwork(i) = 0
               end if
            end select
         end subroutine move_flg

      end subroutine move_extra_info

      end module run_star_extras
      
