module stduse
    use iso_fortran_env, only: real32, real64, input_unit, output_unit, iostat_end
    implicit none

    integer, parameter :: sp     = real32
    integer, parameter :: dp     = real64
    integer, parameter :: stdin  = input_unit
    integer, parameter :: stdout = output_unit
endmodule stduse


module statistics_mod
    use stduse
    implicit none
    private

    public :: vector_avg
    public :: jackknife

    interface vector_avg
        module procedure :: svector_avg
        module procedure :: dvector_avg
        module procedure :: cvector_avg
        module procedure :: zvector_avg
        module procedure :: ivector_avg
    endinterface vector_avg

    interface jackknife
        module procedure :: sjackknife
        module procedure :: djackknife
        module procedure :: cjackknife
        module procedure :: zjackknife
        module procedure :: ijackknife
    endinterface jackknife

    contains

        pure real(sp) function svector_avg(x)
            real(sp), intent(in) :: x(:)
            svector_avg = sum(x) / size(x)
        endfunction svector_avg
        pure real(dp) function dvector_avg(x)
            real(dp), intent(in) :: x(:)
            dvector_avg = sum(x) / size(x)
        endfunction dvector_avg
        pure complex(sp) function cvector_avg(x)
            complex(sp), intent(in) :: x(:)
            cvector_avg = sum(x) / size(x)
        endfunction cvector_avg
        pure complex(dp) function zvector_avg(x) 
            complex(dp), intent(in) :: x(:)
            zvector_avg = sum(x) / size(x)
        endfunction zvector_avg
        pure real(dp) function ivector_avg(x) 
            integer, intent(in) :: x(:)
            ivector_avg = real(sum(x), dp) / size(x)
        endfunction ivector_avg

        pure subroutine sjackknife(x, avg, err)
            real(sp), intent(in)  :: x(:)
            real(sp), intent(out) :: avg, err
            avg = vector_avg(x)
            err = sqrt(sum((x - avg) ** 2) / (size(x)*(size(x)-1)))
        endsubroutine sjackknife
        pure subroutine djackknife(x, avg, err)
            real(dp), intent(in)  :: x(:)
            real(dp), intent(out) :: avg, err
            avg = vector_avg(x)
            err = sqrt(sum((x - avg) ** 2) / (size(x)*(size(x)-1)))
        endsubroutine djackknife
        pure subroutine cjackknife(x, avg, err)
            complex(sp), intent(in)  :: x(:)
            complex(sp), intent(out) :: avg, err
            avg = vector_avg(x)
            err = sqrt(sum((x - avg) ** 2) / (size(x)*(size(x)-1)))
        endsubroutine cjackknife
        pure subroutine zjackknife(x, avg, err)
            complex(dp), intent(in)  :: x(:)
            complex(dp), intent(out) :: avg, err
            avg = vector_avg(x)
            err = sqrt(sum((x - avg) ** 2) / (size(x)*(size(x)-1)))
        endsubroutine zjackknife
        pure subroutine ijackknife(x, avg, err)
            integer , intent(in)  :: x(:)
            real(dp), intent(out) :: avg, err
            avg = vector_avg(x)
            err = sqrt(sum((x - avg) ** 2) / (size(x)*(size(x)-1)))
        endsubroutine ijackknife


        pure function dmatrix_avg(A, dim) result(avg)
            real(dp),           intent(in) :: A(:, :)
            integer , optional, intent(in) :: dim
            integer  :: actualdim
            real(dp) :: avg(size(A, 1))
            if (present(dim)) then
                actualdim = dim
            else
                actualdim = 1
            endif
            avg = sum(A, 2) / size(A, 2)
        endfunction dmatrix_avg


endmodule statistics_mod


module measurementtypes_mod
    use stduse
    use statistics_mod
    implicit none
    private

    public :: scamesr_dp
    public :: scamesr_int
    public :: vecmesr_dp
    !public :: matmesr

    type bincounter
        ! For keeping track of what bin is being filled
        ! and how full it is.
        ! Call self%next() after a measurement for determining where to put the next one
        ! self%setup(nbin, binsize) should immediately be called after initialization
        ! to set up the counter to be used.
        integer :: nbin
        integer :: binsize
        integer :: i   = 1
        integer :: bin = 1
        contains
            procedure :: setup => bincounter_setup
            procedure :: next  => bincounter_next
    endtype bincounter

    type scamesr_dp
        integer :: nbin
        integer :: binsize

        real(dp), allocatable :: bin(:)     ! binsize
        real(dp), allocatable :: binavgs(:) ! nbin
        real(dp) :: avg
        real(dp) :: err

        contains
            procedure :: setup        => scamesr_dp_setup
            procedure :: jackknife    => scamesr_dp_jackknife
            procedure :: avgbin       => scamesr_dp_avgbin
            procedure :: sgnavgadjust => scamesr_dp_sgnavgadjust
            procedure :: unset        => scamesr_dp_unset
    endtype scamesr_dp

    type scamesr_int
        integer :: nbin
        integer :: binsize
        
        integer , allocatable :: bin(:)  ! binsize
        real(dp), allocatable :: binavgs(:) ! nbin
        real(dp) :: avg
        real(dp) :: err

        contains
            procedure :: setup     => scamesr_int_setup
            procedure :: jackknife => scamesr_int_jackknife
            procedure :: avgbin    => scamesr_int_avgbin
            procedure :: unset     => scamesr_int_unset
    endtype scamesr_int

    type vecmesr_dp
        !
        !
        ! bin i = bin(:, i)
        ! Fill column i of bin with vector measurement
        ! bin(j, i) = vector component j in bin i
        !
        integer :: m
        integer :: nbin
        integer :: binsize

        real(dp), allocatable :: bin(:, :)      ! m x binsize
        real(dp), allocatable :: binavgs(:, :)  ! m x nbin
        real(dp), allocatable :: avg(:)         ! m
        real(dp), allocatable :: err(:)         ! m

        contains
            procedure :: setup => vecmesr_dp_setup
    endtype vecmesr_dp

    type matmesr(wp, m, n, nbin, binsize)
        integer, kind :: wp = dp
        integer, len  :: m
        integer, len  :: n
        integer, len  :: nbin
        integer, len  :: binsize

        real(wp) :: bin(m, n, binsize)
        real(wp) :: binavgs(m, n, nbin)
        real(wp) :: avg(m, n)
        real(wp) :: err(m, n)
    endtype matmesr


    contains


        subroutine bincounter_next(self)
            class(bincounter), intent(inout) :: self
            associate(i => self%i, nbin => self%nbin, binsize => self%binsize, bin => self%bin)
                if (i+1 .eq. binsize+1) then
                    i   = 1
                    bin = bin + 1
                else
                    i = i + 1
                endif
            endassociate
        endsubroutine bincounter_next
        subroutine bincounter_setup(self, nbin, binsize)
            class(bincounter), intent(out) :: self
            integer          , intent(in)  :: nbin
            integer          , intent(in)  :: binsize
            self%nbin    = nbin
            self%binsize = binsize
            self%i       = 1
            self%bin     = 1
        endsubroutine bincounter_setup


        subroutine scamesr_dp_setup(self, nbin, binsize)
            class(scamesr_dp), intent(out) :: self
            integer          , intent(in)  :: nbin
            integer          , intent(in)  :: binsize
            self%nbin    = nbin
            self%binsize = binsize
            if (.not. allocated(self%bin)    ) allocate(self%bin(binsize))
            if (.not. allocated(self%binavgs)) allocate(self%binavgs(nbin))
        endsubroutine scamesr_dp_setup
        subroutine scamesr_dp_unset(self)
            class(scamesr_dp), intent(out) :: self
            if (allocated(self%bin)    ) deallocate(self%bin)
            if (allocated(self%binavgs)) deallocate(self%binavgs)
        endsubroutine scamesr_dp_unset
        subroutine scamesr_dp_jackknife(self)
            class(scamesr_dp), intent(inout) :: self
            associate(binavgs => self%binavgs, nbin => self%nbin, avg => self%avg, err => self%err)
                call jackknife(binavgs, nbin, avg, err)
            endassociate
        endsubroutine scamesr_dp_jackknife
        subroutine scamesr_dp_avgbin(self, i)
            ! Does not take into account division by sign average
            ! Later on: binavgs = binavgs / sgnbinavgs
            class(scamesr_dp), intent(inout) :: self
            integer          , intent(in)    :: i
            associate(binavgs => self%binavgs, bin => self%bin, binsize => self%binsize)
                binavgs(i) = vector_avg(bin, binsize)
            endassociate
        endsubroutine scamesr_dp_avgbin
        subroutine scamesr_dp_sgnavgadjust(self, sgnmeas)
            class(scamesr_dp) , intent(inout) :: self
            class(scamesr_int), intent(in)    :: sgnmeas
            associate(binavgs => self%binavgs, sgnbinavgs => sgnmeas%binavgs)
                binavgs = binavgs / sgnbinavgs
            endassociate
        endsubroutine scamesr_dp_sgnavgadjust


        subroutine scamesr_int_setup(self, nbin, binsize)
            class(scamesr_int), intent(out) :: self
            integer           , intent(in)  :: nbin
            integer           , intent(in)  :: binsize
            self%nbin    = nbin
            self%binsize = binsize
            if (.not. allocated(self%bin)    ) allocate(self%bin(binsize))
            if (.not. allocated(self%binavgs)) allocate(self%binavgs(nbin))
        endsubroutine scamesr_int_setup
        subroutine scamesr_int_unset(self)
            class(scamesr_int), intent(out) :: self
            if (allocated(self%bin)    ) deallocate(self%bin)
            if (allocated(self%binavgs)) deallocate(self%binavgs)
        endsubroutine scamesr_int_unset
        subroutine scamesr_int_jackknife(self)
            class(scamesr_int), intent(inout) :: self
            associate(binavgs => self%binavgs, nbin => self%nbin, avg => self%avg, err => self%err)
                call jackknife(binavgs, nbin, avg, err)
            endassociate
        endsubroutine scamesr_int_jackknife
        subroutine scamesr_int_avgbin(self, i)
            class(scamesr_int), intent(inout) :: self
            integer           , intent(in)    :: i
            associate(binavgs => self%binavgs, bin => self%bin, binsize => self%binsize)
                binavgs(i) = vector_avg(bin, binsize)
            endassociate
        endsubroutine scamesr_int_avgbin


        subroutine vecmesr_dp_setup(self, m, nbin, binsize)
            class(vecmesr_dp), intent(out) :: self
            integer          , intent(in)  :: m
            integer          , intent(in)  :: nbin
            integer          , intent(in)  :: binsize
            if (.not. allocated(self%bin)    ) allocate(self%bin(m, binsize))
            if (.not. allocated(self%binavgs)) allocate(self%binavgs(m, nbin))
            if (.not. allocated(self%avg)    ) allocate(self%avg(m))
            if (.not. allocated(self%err)    ) allocate(self%err(m))
        endsubroutine vecmesr_dp_setup
        subroutine vecmesr_dp_unset(self)
            class(vecmesr_dp), intent(out) :: self
            if (allocated(self%bin)    ) deallocate(self%bin)
            if (allocated(self%binavgs)) deallocate(self%binavgs)
            if (allocated(self%avg)    ) deallocate(self%avg)
            if (allocated(self%err)    ) deallocate(self%err)
        endsubroutine vecmesr_dp_unset
        subroutine vecmesr_dp_avgbin(self, i)
            ! Does not take into account division by sign average
            ! Later on: binavgs = binavgs / sgnbinavgs
            class(scamesr_dp), intent(inout) :: self
            integer          , intent(in)    :: i
            associate(binavgs => self%binavgs, bin => self%bin, binsize => self%binsize)
                binavgs(i) = vector_avg(bin, binsize)
            endassociate
        endsubroutine vecmesr_dp_avgbin

        
endmodule measurementtypes_mod