program hopping
    use, intrinsic :: ISO_Fortran_env
    implicit none

    integer, parameter :: msize = 1000000
    integer, parameter :: kindr = selected_real_kind(15, 307)
    integer, allocatable :: nidv1neigh(:,:), nidv2neigh(:,:)
    integer, allocatable :: idv1delete(:), idv2delete(:)
    character(100) :: infile1, outfile1, v1neigh, v2neigh
    character(100) :: tddata, v1data, v2data, ohdata, switched
    character(2), allocatable :: atom(:), atomtd(:), atomv1(:), atomv2(:), atomoh(:)
    character(2), allocatable :: atomv1neigh(:), atomv2neigh(:)
    real(kindr), allocatable :: x(:), y(:), z(:)
    real(kindr), allocatable :: xtd(:), ytd(:), ztd(:)
    real(kindr), allocatable :: xv1(:), yv1(:), zv1(:)
    real(kindr), allocatable :: xv2(:), yv2(:), zv2(:)
    real(kindr), allocatable :: xoh(:), yoh(:), zoh(:)
    real(kindr), allocatable :: xv1neigh(:), yv1neigh(:), zv1neigh(:)
    real(kindr), allocatable :: xv2neigh(:), yv2neigh(:), zv2neigh(:)
    real(kindr), allocatable :: v1_x(:,:), v1_y(:,:), v1_z(:,:)
    real(kindr), allocatable :: v2_x(:,:), v2_y(:,:), v2_z(:,:)

    real(kindr) :: xx, yy, zz, tempx, tempy, tempz, d
    real(kindr) :: k_b, temperature, hop_probability, r, u
    real(kindr) :: total_p, cumulative_prob
    integer :: ios, na, ntd, nv1, nv2, noh, nv1neigh, nv2neigh
    integer :: ntot, numattempts, nrand, i, j, k, m, ii, jj
    integer :: total_paths, nv1v2, selected_path

    ! 10 activation energies for distinct migration pathways
    real(kindr), dimension(10) :: activation_energies = (/ &
        0.79, 0.79, 0.76, 0.83, &       ! V1 pathways (I, IV, II, III)
        0.76, 0.76, 0.76, 0.79, 0.79, 0.79 /)  ! V2 pathways (V-X)

    real(kindr), dimension(10) :: probabilities, relative_fractions, accumulated_probabilities
    integer, dimension(10) :: path_counts

    path_counts = 0
    accumulated_probabilities = 0.0

    k_b = 8.617333262145e-5  ! Boltzmann constant in eV/K
    temperature = 300.0      ! Temperature in Kelvin

    infile1 = '100_100superblock.xyz'
    outfile1 = 'outfile1.data'
    v1neigh = 'v1neigh.xyz'
    v2neigh = 'v2neigh.xyz'
    tddata = 'tddata.xyz'
    v1data = 'v1data.xyz'
    v2data = 'v2data.xyz'
    ohdata = 'ohdata.xyz'
    switched = 'switched.xyz'

    na = 680002

    ! Allocate arrays
    allocate(nidv1neigh(msize, 6), nidv2neigh(msize, 6))
    allocate(idv1delete(msize), idv2delete(msize))
    allocate(atom(msize), atomtd(msize), atomv1(msize), atomv2(msize), atomoh(msize))
    allocate(atomv1neigh(msize), atomv2neigh(msize))
    allocate(x(msize), y(msize), z(msize))
    allocate(xtd(msize), ytd(msize), ztd(msize))
    allocate(xv1(msize), yv1(msize), zv1(msize))
    allocate(xv2(msize), yv2(msize), zv2(msize))
    allocate(xoh(msize), yoh(msize), zoh(msize))
    allocate(xv1neigh(msize), yv1neigh(msize), zv1neigh(msize))
    allocate(xv2neigh(msize), yv2neigh(msize), zv2neigh(msize))
    allocate(v1_x(msize, 6), v1_y(msize, 6), v1_z(msize, 6))
    allocate(v2_x(msize, 6), v2_y(msize, 6), v2_z(msize, 6))

    ! Read input file
    open(unit = 10, file = infile1, status = 'old', iostat = ios)
    if (ios /= 0) then
        print *, "Error opening file: ", infile1
        stop
    endif
    read(10, *)
    read(10, *)
    do i = 1, na - 2
        read(10, *) atom(i), x(i), y(i), z(i)
    end do
    close(10)

    open(unit = 21, file = outfile1, status = 'replace')
    do i = 1, na - 2
        write(21, *) atom(i), x(i), y(i), z(i)
    end do
    close(21)

    ! Count atom types
    ntd = 0; nv1 = 0; nv2 = 0; noh = 0
    do i = 1, na - 1
        if (atom(i) .EQ. 'Mn') then
            ntd = ntd + 1
            atomtd(ntd) = atom(i)
            xtd(ntd) = x(i); ytd(ntd) = y(i); ztd(ntd) = z(i)
        else if (atom(i) .EQ. 'B ') then
            nv1 = nv1 + 1
            atomv1(nv1) = atom(i)
            xv1(nv1) = x(i); yv1(nv1) = y(i); zv1(nv1) = z(i)
        else if (atom(i) .EQ. 'C ') then
            nv2 = nv2 + 1
            atomv2(nv2) = atom(i)
            xv2(nv2) = x(i); yv2(nv2) = y(i); zv2(nv2) = z(i)
        else if (atom(i) .EQ. 'Fe') then
            noh = noh + 1
            atomoh(noh) = atom(i)
            xoh(noh) = x(i); yoh(noh) = y(i); zoh(noh) = z(i)
        end if
    end do
    write(*, *) 'Number of Td Fe atoms = ', ntd
    write(*, *) 'Number of V1 vacancies = ', nv1
    write(*, *) 'Number of V2 vacancies = ', nv2
    write(*, *) 'Number of Oh Fe atoms = ', noh

    ! V1 neighbor search
    nv1neigh = 0
    do i = 1, nv1
        k = 0
        do j = 1, noh
            xx = abs(xv1(i) - xoh(j))
            yy = abs(yv1(i) - yoh(j))
            zz = abs(zv1(i) - zoh(j))
            d = sqrt(xx * xx + yy * yy + zz * zz)
            if (d .LE. 3.5) then
                nv1neigh = nv1neigh + 1
                k = k + 1
                atomv1neigh(nv1neigh) = atomoh(j)
                xv1neigh(nv1neigh) = xoh(j)
                yv1neigh(nv1neigh) = yoh(j)
                zv1neigh(nv1neigh) = zoh(j)
                v1_x(i, k) = xv1neigh(nv1neigh)
                v1_y(i, k) = yv1neigh(nv1neigh)
                v1_z(i, k) = zv1neigh(nv1neigh)
                nidv1neigh(i, k) = j
            end if
        end do
    end do
    write(*, *) 'Number of V1 neighbors = ', nv1neigh

    ntot = nv1 + nv1neigh
    open(unit = 30, file = v1neigh, status = 'replace')
    write(30, *) ntot
    write(30, *) 'V1 with neighbors'
    do j = 1, nv1neigh
        write(30, *) atomv1neigh(j), xv1neigh(j), yv1neigh(j), zv1neigh(j)
    end do
    do i = 1, nv1
        write(30, *) atomv1(i), xv1(i), yv1(i), zv1(i)
    end do
    close(30)

    ! V2 neighbor search
    nv2neigh = 0
    do i = 1, nv2
        k = 0
        do j = 1, noh
            xx = abs(xv2(i) - xoh(j))
            yy = abs(yv2(i) - yoh(j))
            zz = abs(zv2(i) - zoh(j))
            d = sqrt(xx * xx + yy * yy + zz * zz)
            if (d .LE. 3.5) then
                nv2neigh = nv2neigh + 1
                k = k + 1
                atomv2neigh(nv2neigh) = atomoh(j)
                xv2neigh(nv2neigh) = xoh(j)
                yv2neigh(nv2neigh) = yoh(j)
                zv2neigh(nv2neigh) = zoh(j)
                v2_x(i, k) = xv2neigh(nv2neigh)
                v2_y(i, k) = yv2neigh(nv2neigh)
                v2_z(i, k) = zv2neigh(nv2neigh)
                nidv2neigh(i, k) = j
            end if
        end do
    end do
    write(*, *) 'Number of V2 neighbors = ', nv2neigh

    ntot = nv2 + nv2neigh
    open(unit = 30, file = v2neigh, status = 'replace')
    write(30, *) ntot
    write(30, *) 'V2 with neighbors'
    do j = 1, nv2neigh
        write(30, *) atomv2neigh(j), xv2neigh(j), yv2neigh(j), zv2neigh(j)
    end do
    do i = 1, nv2
        write(30, *) atomv2(i), xv2(i), yv2(i), zv2(i)
    end do
    close(30)

    ! KMC hopping loop with multi-pathway selection
    nv1v2 = nv1 + nv2
    numattempts = int(nv1v2 * 10)
    do k = 1, numattempts
        call random_number(u)
        nrand = 1 + INT((nv1v2) * u)

        if (nrand .LE. nv1) then
            call random_number(u)
            i = 1 + FLOOR((nv1) * u)
            if ((v1_x(i, 4) .NE. 0) .AND. any(idv1delete .NE. i)) then
                call random_number(u)
                j = 1 + FLOOR((4) * u)
                m = nidv1neigh(i, j)
                do ii = 1, 4
                    probabilities(ii) = exp(-activation_energies(ii) / (k_b * temperature))
                end do
                probabilities(5:10) = 0.0
            else
                cycle
            end if
        else
            call random_number(u)
            i = 1 + FLOOR((nv2) * u)
            if ((v2_x(i, 6) .NE. 0) .AND. any(idv2delete .NE. i)) then
                call random_number(u)
                j = 1 + FLOOR((6) * u)
                m = nidv2neigh(i, j)
                probabilities(1:4) = 0.0
                do ii = 5, 10
                    probabilities(ii) = exp(-activation_energies(ii) / (k_b * temperature))
                end do
            else
                cycle
            end if
        end if

        ! Normalize and select path
        total_p = sum(probabilities)
        probabilities = probabilities / total_p
        accumulated_probabilities = accumulated_probabilities + probabilities

        call random_number(u)
        cumulative_prob = 0.0
        do ii = 1, 10
            cumulative_prob = cumulative_prob + probabilities(ii)
            if (u < cumulative_prob) then
                selected_path = ii
                exit
            end if
        end do

        path_counts(selected_path) = path_counts(selected_path) + 1
        hop_probability = probabilities(selected_path)

        call random_number(r)
        if (r < hop_probability) then
            if (nrand .LE. nv1) then
                tempx = xv1(i); tempy = yv1(i); tempz = zv1(i)
                xv1(i) = v1_x(i, j); yv1(i) = v1_y(i, j); zv1(i) = v1_z(i, j)
                xoh(m) = tempx; yoh(m) = tempy; zoh(m) = tempz
                idv1delete(k) = i
            else
                tempx = xv2(i); tempy = yv2(i); tempz = zv2(i)
                xv2(i) = v2_x(i, j); yv2(i) = v2_y(i, j); zv2(i) = v2_z(i, j)
                xoh(m) = tempx; yoh(m) = tempy; zoh(m) = tempz
                idv2delete(k) = i
            end if

            do ii = 1, nv1
                do jj = 1, 6
                    if (nidv1neigh(ii, jj) .EQ. m) nidv1neigh(ii, jj) = 0
                end do
            end do
            do ii = 1, nv2
                do jj = 1, 6
                    if (nidv2neigh(ii, jj) .EQ. m) nidv2neigh(ii, jj) = 0
                end do
            end do
        end if
    end do

    ! Normalize and report pathway statistics
    total_p = sum(accumulated_probabilities)
    accumulated_probabilities = accumulated_probabilities / total_p
    total_paths = sum(path_counts)
    do i = 1, 10
        relative_fractions(i) = real(path_counts(i), kindr) / real(total_paths, kindr)
    end do

    write(*, *) 'Pathway Results:'
    write(*, *) '------------------------------------------------------------'
    write(*, *) 'Path | Ea (eV) | Rel. Probability | Count | Rel. Fraction'
    write(*, *) '------------------------------------------------------------'
    do i = 1, 10
        write(*, '(I4, " | ", F7.2, " | ", F16.6, " | ", I5, " | ", F13.6)') &
            i, activation_energies(i), accumulated_probabilities(i), &
            path_counts(i), relative_fractions(i)
    end do
    write(*, *) '------------------------------------------------------------'

    ! Write final structure
    ntot = nv1 + nv2 + noh
    open(unit = 30, file = switched, status = 'replace')
    write(30, *) ntot
    write(30, *) 'V1 and V2 with Oh Fe'
    do j = 1, noh
        write(30, *) 'Fe', xoh(j), yoh(j), zoh(j)
    end do
    do i = 1, nv1
        write(30, *) atomv1(i), xv1(i), yv1(i), zv1(i)
    end do
    do i = 1, nv2
        write(30, *) atomv2(i), xv2(i), yv2(i), zv2(i)
    end do
    close(30)

    write(*, *) 'number of v1 switches = ', maxval(idv1delete)
    write(*, *) 'number of v2 switches = ', maxval(idv2delete)

    ! Deallocate
    deallocate(nidv1neigh, nidv2neigh, idv1delete, idv2delete)
    deallocate(atom, atomtd, atomv1, atomv2, atomoh)
    deallocate(atomv1neigh, atomv2neigh)
    deallocate(x, y, z, xtd, ytd, ztd)
    deallocate(xv1, yv1, zv1, xv2, yv2, zv2, xoh, yoh, zoh)
    deallocate(xv1neigh, yv1neigh, zv1neigh, xv2neigh, yv2neigh, zv2neigh)
    deallocate(v1_x, v1_y, v1_z, v2_x, v2_y, v2_z)

end program hopping
