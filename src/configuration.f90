module configuration
    use kinds, only: dp
    implicit none
    private

    public :: toolkit_config
    public :: load_config_file
    public :: parse_command_line
    public :: print_help
    public :: print_config

    type :: toolkit_config
        character(len=256) :: input_file = ''
        character(len=256) :: saddles_file = 'saddles.dat'
        character(len=256) :: separatrices_file = 'separatrices.dat'
        character(len=256) :: diagnostics_file = 'diagnostics.txt'

        integer :: layer_id = 1

        ! Negative values mean "automatic".
        real(dp) :: root_tol = -1.0_dp
        real(dp) :: merge_tol = -1.0_dp
        real(dp) :: level_tol = -1.0_dp
        real(dp) :: radius_factor = 2.0_dp
        real(dp) :: connect_tol = -1.0_dp
        real(dp) :: min_cos = 0.90_dp

        logical :: write_diagnostics = .true.
        logical :: verbose = .true.
    end type toolkit_config

contains

    subroutine load_config_file(filename, cfg)
        character(len=*), intent(in) :: filename
        type(toolkit_config), intent(inout) :: cfg

        integer :: unit, ios
        character(len=256) :: input_file, saddles_file, separatrices_file, diagnostics_file
        integer :: layer_id
        real(dp) :: root_tol, merge_tol, level_tol, radius_factor, connect_tol, min_cos
        logical :: write_diagnostics, verbose

        namelist /separatrix/ input_file, saddles_file, separatrices_file, diagnostics_file, &
                              layer_id, root_tol, merge_tol, level_tol, radius_factor, &
                              connect_tol, min_cos, write_diagnostics, verbose

        input_file = cfg%input_file
        saddles_file = cfg%saddles_file
        separatrices_file = cfg%separatrices_file
        diagnostics_file = cfg%diagnostics_file
        layer_id = cfg%layer_id
        root_tol = cfg%root_tol
        merge_tol = cfg%merge_tol
        level_tol = cfg%level_tol
        radius_factor = cfg%radius_factor
        connect_tol = cfg%connect_tol
        min_cos = cfg%min_cos
        write_diagnostics = cfg%write_diagnostics
        verbose = cfg%verbose

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,'(A,A)') 'ERROR: cannot open config file: ', trim(filename)
            error stop 2
        end if

        read(unit, nml=separatrix, iostat=ios)
        close(unit)

        if (ios /= 0) then
            write(*,'(A,A)') 'ERROR: cannot parse config file: ', trim(filename)
            error stop 2
        end if

        cfg%input_file = input_file
        cfg%saddles_file = saddles_file
        cfg%separatrices_file = separatrices_file
        cfg%diagnostics_file = diagnostics_file
        cfg%layer_id = layer_id
        cfg%root_tol = root_tol
        cfg%merge_tol = merge_tol
        cfg%level_tol = level_tol
        cfg%radius_factor = radius_factor
        cfg%connect_tol = connect_tol
        cfg%min_cos = min_cos
        cfg%write_diagnostics = write_diagnostics
        cfg%verbose = verbose
    end subroutine load_config_file


    subroutine parse_command_line(cfg)
        type(toolkit_config), intent(inout) :: cfg

        integer :: argc, i
        character(len=256) :: arg, value, config_file
        logical :: has_config

        argc = command_argument_count()
        has_config = .false.
        config_file = ''

        ! Pass 1: find config file so command-line values can override it
        ! regardless of argument order.
        i = 1
        do while (i <= argc)
            call get_command_argument(i, arg)
            if (trim(arg) == '--config') then
                if (i == argc) error stop '--config requires a filename'
                call get_command_argument(i+1, config_file)
                has_config = .true.
                i = i + 2
            else
                i = i + 1
            end if
        end do

        if (has_config) call load_config_file(trim(config_file), cfg)

        ! Pass 2: apply CLI overrides.
        i = 1
        do while (i <= argc)
            call get_command_argument(i, arg)

            select case (trim(arg))
            case ('--config')
                i = i + 2

            case ('--input')
                call require_value(i, argc, value)
                cfg%input_file = value
                i = i + 2

            case ('--layer')
                call require_value(i, argc, value)
                read(value,*) cfg%layer_id
                i = i + 2

            case ('--root-tol')
                call require_value(i, argc, value)
                read(value,*) cfg%root_tol
                i = i + 2

            case ('--merge-tol')
                call require_value(i, argc, value)
                read(value,*) cfg%merge_tol
                i = i + 2

            case ('--level-tol')
                call require_value(i, argc, value)
                read(value,*) cfg%level_tol
                i = i + 2

            case ('--radius-factor')
                call require_value(i, argc, value)
                read(value,*) cfg%radius_factor
                i = i + 2

            case ('--connect-tol')
                call require_value(i, argc, value)
                read(value,*) cfg%connect_tol
                i = i + 2

            case ('--min-cos')
                call require_value(i, argc, value)
                read(value,*) cfg%min_cos
                i = i + 2

            case ('--saddles-output')
                call require_value(i, argc, value)
                cfg%saddles_file = value
                i = i + 2

            case ('--separatrices-output')
                call require_value(i, argc, value)
                cfg%separatrices_file = value
                i = i + 2

            case ('--diagnostics-output')
                call require_value(i, argc, value)
                cfg%diagnostics_file = value
                i = i + 2

            case ('--no-diagnostics')
                cfg%write_diagnostics = .false.
                i = i + 1

            case ('--quiet')
                cfg%verbose = .false.
                i = i + 1

            case ('--help', '-h')
                call print_help()
                stop

            case default
                if (len_trim(arg) > 0 .and. arg(1:1) == '-') then
                    write(*,'(A,A)') 'ERROR: unknown option: ', trim(arg)
                    call print_help()
                    error stop 2
                end if

                ! Backward-compatible positional input filename.
                if (len_trim(cfg%input_file) == 0) then
                    cfg%input_file = arg
                else
                    write(*,'(A,A)') 'ERROR: unexpected positional argument: ', trim(arg)
                    error stop 2
                end if
                i = i + 1
            end select
        end do

        if (len_trim(cfg%input_file) == 0) then
            call print_help()
            error stop 'No input field was specified'
        end if

        if (cfg%radius_factor <= 0.0_dp) error stop 'radius_factor must be positive'
        if (cfg%min_cos <= 0.0_dp .or. cfg%min_cos > 1.0_dp) &
            error stop 'min_cos must be in (0,1]'
    end subroutine parse_command_line


    subroutine require_value(i, argc, value)
        integer, intent(in) :: i, argc
        character(len=256), intent(out) :: value

        if (i >= argc) error stop 'Command-line option requires a value'
        call get_command_argument(i+1, value)
    end subroutine require_value


    subroutine print_config(cfg)
        type(toolkit_config), intent(in) :: cfg

        write(*,'(A)') '--- Separatrix Toolkit configuration ---'
        write(*,'(A,A)') 'input_file          = ', trim(cfg%input_file)
        write(*,'(A,I0)') 'layer_id            = ', cfg%layer_id
        write(*,'(A,ES12.4)') 'root_tol            = ', cfg%root_tol
        write(*,'(A,ES12.4)') 'merge_tol           = ', cfg%merge_tol
        write(*,'(A,ES12.4)') 'level_tol           = ', cfg%level_tol
        write(*,'(A,F8.3)') 'radius_factor       = ', cfg%radius_factor
        write(*,'(A,ES12.4)') 'connect_tol         = ', cfg%connect_tol
        write(*,'(A,F8.3)') 'min_cos             = ', cfg%min_cos
        write(*,'(A,A)') 'saddles_file        = ', trim(cfg%saddles_file)
        write(*,'(A,A)') 'separatrices_file   = ', trim(cfg%separatrices_file)
        write(*,'(A,A)') 'diagnostics_file    = ', trim(cfg%diagnostics_file)
        write(*,'(A)') '----------------------------------------'
    end subroutine print_config


    subroutine print_help()
        write(*,'(A)') 'Separatrix Toolkit'
        write(*,'(A)') ''
        write(*,'(A)') 'Usage:'
        write(*,'(A)') '  separatrix_toolkit FIELD.dat [options]'
        write(*,'(A)') '  separatrix_toolkit --config config.nml [options]'
        write(*,'(A)') ''
        write(*,'(A)') 'Options:'
        write(*,'(A)') '  --input FILE'
        write(*,'(A)') '  --config FILE'
        write(*,'(A)') '  --layer N'
        write(*,'(A)') '  --root-tol X'
        write(*,'(A)') '  --merge-tol X'
        write(*,'(A)') '  --level-tol X'
        write(*,'(A)') '  --radius-factor X'
        write(*,'(A)') '  --connect-tol X'
        write(*,'(A)') '  --min-cos X'
        write(*,'(A)') '  --saddles-output FILE'
        write(*,'(A)') '  --separatrices-output FILE'
        write(*,'(A)') '  --diagnostics-output FILE'
        write(*,'(A)') '  --no-diagnostics'
        write(*,'(A)') '  --quiet'
        write(*,'(A)') '  --help'
        write(*,'(A)') ''
        write(*,'(A)') 'Negative tolerance values select automatic defaults.'
    end subroutine print_help

end module configuration
