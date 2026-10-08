program separatrix_toolkit
    use separatrix_types, only: grid_field, saddle_point, separatrix_network
    use grid_io, only: read_xyz_grid
    use critical_points, only: find_all_saddles
    use saddle_refinement, only: refine_saddles_quadratic
    use separatrix_builder, only: build_all_separatrices
    use output_writer, only: write_saddles, write_separatrices
    use configuration, only: toolkit_config, parse_command_line, print_config
    use diagnostics, only: write_diagnostics_report, print_topology_summary
    implicit none

    type(toolkit_config) :: cfg
    type(grid_field) :: field
    type(saddle_point), allocatable :: saddles(:)
    type(separatrix_network), allocatable :: networks(:)

    integer :: nsaddles, i

    call parse_command_line(cfg)

    if (cfg%verbose) call print_config(cfg)

    call read_xyz_grid(trim(cfg%input_file), field)

    call find_all_saddles(field, saddles, nsaddles, &
                          root_tol=cfg%root_tol, merge_tol=cfg%merge_tol)

    call refine_saddles_quadratic(field, saddles)

    nsaddles = size(saddles)
    do i = 1, nsaddles
        saddles(i)%id = i
    end do

    call build_all_separatrices(field, saddles, networks, &
                                level_tol=cfg%level_tol, &
                                radius_factor=cfg%radius_factor, &
                                connect_tol=cfg%connect_tol, &
                                min_cos=cfg%min_cos)

    call write_saddles(trim(cfg%saddles_file), cfg%layer_id, saddles)
    call write_separatrices(trim(cfg%separatrices_file), cfg%layer_id, networks)

    if (cfg%write_diagnostics) then
        call write_diagnostics_report(trim(cfg%diagnostics_file), field, saddles, networks)
    end if

    if (cfg%verbose) then
        call print_topology_summary(saddles, networks)
        write(*,'(A,A)') 'Saddles: ', trim(cfg%saddles_file)
        write(*,'(A,A)') 'Separatrices: ', trim(cfg%separatrices_file)
        if (cfg%write_diagnostics) write(*,'(A,A)') 'Diagnostics: ', trim(cfg%diagnostics_file)
    end if

end program separatrix_toolkit
