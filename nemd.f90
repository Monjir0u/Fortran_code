! 変数名のルール
! num → 数字, 数 (number)
! mol → 分子 (molecule)
! temp → 温度 (temepature)
! time → 時間 (time)
! step → ステップ (step)
! up → 上壁面 (upper layer)
! lw → 下壁面 (lower layer)
! coef → 係数 (coefficient)
! heat → 熱量 (heat)
! pot → ポテンシャルエネルギー (potential energy)
! kin → 運動エネルギー (kinetic energy)
! len → 長さ (length)

! 定数は大文字で定義

module parameters
    implicit none
    !--------------- よく変更する---------------------!
    ! 計算時間
    double precision, parameter :: time_scaling = 0.1d0 ! スケーリング時間 [ns]
    double precision, parameter :: time_relax   = 1.0d0 ! 緩和計算時間 [ns]
    double precision, parameter :: time_measure = 10.0d0 ! データ計測時間 [ns]
    double precision, parameter :: time_all = time_scaling + time_relax + time_measure ! 全計算時間
    
    double precision, parameter :: dt = 1.00d0 ! 無次元時間ステップ(無次元，有次元の場合fs)
    double precision, parameter :: tau = 1.00d0 ! 測定間隔 [fs]

    double precision, parameter :: TEMP_Ar = 100d0 ! 系内（目標）温度  [K]
    double precision, parameter :: INTER_STG = 0.1d0 ! 相互作用強さ (interaction strength)

    character(len=20) :: dir_name = 'a0.1'  ! ファイル名
    !------------------------------------------------!

    ! ステップ
    integer, parameter :: step_scaling = int(time_scaling/dt*1.0d+6) ! スケーリングステップ
    integer, parameter :: step_relax   = int(time_relax  /dt*1.0d+6) ! 緩和計算ステップ
    integer, parameter :: step_measure = int(time_measure/dt*1.0d+6) ! データ計測ステップ
    integer, parameter :: step_all = step_scaling + step_relax + step_measure ! 最大ステップ数

    ! 分子の数，配置
    integer, parameter :: TYPMOL = 3 ! 分子の種類数
    integer, parameter :: NUM_X(TYPMOL) = [20,10,20]
    integer, parameter :: NUM_Y(TYPMOL) = [10, 8,10]
    integer, parameter :: NUM_Z(TYPMOL) = [ 6,25, 6]
    integer, parameter :: NUM_MOL(TYPMOL) = [NUM_X(1)*NUM_Y(1)*NUM_Z(1), NUM_X(2)*NUM_Y(2)*NUM_Z(2), NUM_X(3)*NUM_Y(3)*NUM_Z(3)] ! 各分子の数
    integer, parameter :: NUM_MOL_ALL = NUM_MOL(1) + NUM_MOL(2) + NUM_MOL(3)
    
    ! 境界条件
    double precision, parameter :: ST_DIST(TYPMOL) = [3.92d0, 4.8d0, 3.92d0] ! 格子定数(無次元)[Å] (stndard distance)
    double precision, parameter :: BND_LEN_X0 = ST_DIST(1) * NUM_Y(1) ! x方向の周期境界長さ(無次元) (boundary length)
    double precision, parameter :: BND_LEN_Y0 = ST_DIST(1) * NUM_Y(1) ! y方向の周期境界長さ(無次元）
    double precision, parameter :: BND_LEN_Z0 = 87.0d0 ! z方向の周期境界長さ(無次元）
    double precision, parameter :: BND_LEN0(TYPMOL) = [BND_LEN_X0, BND_LEN_Y0, BND_LEN_Z0]
    double precision, parameter :: THICK_Pt_UP = ST_DIST(1)*(NUM_Z(1)*0.5d0+0.25d0)
    double precision, parameter :: THICK_Pt_LW = ST_DIST(3)*(NUM_Z(3)*0.5d0+0.25d0)
    double precision, parameter :: THICK_Ar = BND_LEN_Z0 - THICK_Pt_UP - THICK_Pt_LW
    double precision, parameter :: THICK(TYPMOL) = [THICK_Pt_UP, THICK_Ar, THICK_Pt_LW]

    ! 定数
    double precision, parameter :: CUTOFF = 3.300d0 ! カットオフ長さ/σ
    double precision, parameter :: AVOGA = 6.022d+23 ! アボガドロ数
    double precision, parameter :: BOLTZ = 1.3806662d-23 ! ボルツマン定数 [J/K]
    double precision, parameter :: PI = 3.141592654d0 ! 円周率

    ! 分子の質量
    !double precision, parameter :: bunsi(TYPMOL) = [195.084d-3, 39.950d-3, 195.084d-3] ! 分子の質量  kg/mol   
    double precision, parameter :: MASS(TYPMOL) = [32.395d0, 6.6340d0, 32.395d0] ! 分子の質量（無次元） * 10d-26 [kg/個]
    
    ! Lennard-Jonesパラメータ
    double precision, parameter :: SIG(4) = [2.540d0, 3.400d0, 2.540d0, 2.970d0]  ! σ(無次元) *1.0d-10
    double precision, parameter :: EPS(4) = [109.2d-5, 1.666d-5, 109.2d-5, 13.49d-5] ! ε(無次元) *1.0d-16

    ! 粒子登録法
    integer, parameter :: step_update = 40 ! 更新ステップ

    ! Langevin法
    double precision, parameter :: TEMP_Langevin(TYPMOL) = [130d0, 0d0, 70d0] ! Langevin法を用いるPtの温度  真ん中は使わない
    double precision, parameter :: DIRAC = 1.054571817d-34 ! ディラック定数 [J･s]
    double precision, parameter :: DEBTMP = 240d0 ! Debye温度 [K]
    double precision, parameter :: OMEGA = 3.14212728482d+13 ! BOLTZ * DEBTMP / DIRAC * 1.0d-11 ! Debye定数 (有次元)
    double precision, parameter :: DAMP =  5.32967075080d-12 ! MASS(1) * PI * OMEGA / 6.000d0 ! ダンパーの減衰係数 (有次元)
    
    ! 熱流束
    double precision, parameter :: AREA_Pt = ST_DIST(1)*ST_DIST(1)*int(NUM_X(1)*0.5)*NUM_Y(1)
    
    ! 温度分布
    integer, parameter :: NUM_DIV_Ar = 15 ! Arの温度分布の分割数
    double precision :: LEN_DIV_Ar = THICK(2)/NUM_DIV_Ar ! Arの温度分布の分割距離

end module parameters

! 変数
module variable
    use parameters
    implicit none
    integer :: step_now ! 現在のステップ数

    ! 分子間力
    double precision :: coef_force(4) = 0.0d0 ! 力の計算の係数

    ! カットオフ
    double precision :: cutoff(3) ! ポテンシャルのカットオフ長さx,y,z方向
    double precision :: bnd_len(3) ! ポテンシャルのカットオフ長さx,y,z方向，x,y,z方向の周期長さ

    ! Langevin法
    double precision :: force_rand(NUM_MOL(1),TYPMOL,3)   ! ランダム力用
    double precision :: force_damp(NUM_MOL(1),TYPMOL,3)   ! ダンパー力用
    double precision :: rnd_2
    logical :: isOdd = .true. ! 乱数のsinとcosを交互に出すためのフラグ

    ! 熱流束
    double precision :: force_inter(NUM_MOL(1),TYPMOL,3) ! 熱流束を計算するための相互作用力 z方向のみを使う
    double precision :: heat_phantom(TYPMOL)! = 0.0d0 ! Phantom層からの熱輸送量
    double precision :: heat_interface(TYPMOL)! = 0.0d0 ! 固液界面での熱輸送量

    double precision :: temp_layer_Pt(NUM_Z(1),TYPMOL) ! Ptの層ごとの温度
    double precision :: temp_layer_Ar(NUM_DIV_Ar) ! z方向に分割した領域内のArの温度

    contains

    !----------- 関数定義
    ! 乱数生成用の関数
    function Random() result(rnd_) ! 呼び出す時に何度も計算しなくていいように下の関数と分けた
        use parameters
        implicit none
        double precision :: rnd_, u1, u2
        
        ! 呼び出される回数が偶数か奇数かによってsinかcosかを使い分ける
        if(isOdd) then
            call random_number(u1)
            call random_number(u2)
            rnd_2 = sqrt(-2.0d0 * log(u1)) * sin(2.0d0 * PI * u2)
            rnd_  = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * PI * u2)
            isOdd = .false.
        else
           rnd_ = rnd_2
           isOdd = .true.
        end if

    end function random
    
    ! 設定温度を引数とし、ランダム力計算用の標準偏差を出力する関数
    function getStddev(T_) result(stddev_)
        use parameters
        implicit none
        double precision, intent(in) :: T_
        double precision :: stddev_

        stddev_ = dsqrt(2.000d0 * DAMP * BOLTZ * T_ / dt * 1.0d+33)  ! 無次元 stddev -> 有次元では10^9　この辺の無次元化はかなり複雑になっているので、あまりいじらない方がいいです

    end function getStddev
end module variable

! 触りたての最初の方に、同じ処理を省けるように構造体を作りました
! イメージとしては、分子の種類（上Pt, Ar, 下Pt）× 分子の番号の二次元配列があって、その要素一つ一つが位置や速度などの情報を持っている感じ
! ですがコードが長くなるのであまりおすすめしません
! おすすめは分子の種類 × 分子の番号の二次元配列を位置や速度それぞれで作成する方法です　　pos(kind, number)みたいな

module molecules_struct
    use parameters
    implicit none

    ! 構造体の定義
    type :: mol_info
        double precision :: pos(3)
        double precision :: vel(3), vtmp(3)
        double precision :: acc(3)
        double precision :: pot, kin
    end type mol_info

    type mol_typ
        type(mol_info), allocatable :: mol(:) ! 分子の番号を区別
    end type mol_typ

    type(mol_typ), dimension(TYPMOL) :: typ  ! 分子の種類を区別　上Pt, 中Ar, 下Pt

end module molecules_struct

! PVWINという可視化ソフト用
! module forPVwin
!     use parameters
!     implicit none
!     integer, parameter :: NUM_MOL_ALL = NUM_MOL(1) + NUM_MOL(2) + NUM_MOL(3)
!     integer, parameter :: moltype = 1
!     integer, parameter :: ndat = int(step_all/100)
!     integer, parameter :: ntime0 = 0
!     integer, parameter :: ndt = 1
! end module forPVwin
! PVWINはおすすめしません

program main
    use parameters
    use variable
    use molecules_struct
    ! use forPVwin
    implicit none
    double precision :: time
    integer :: i, j

    character(len=100) :: filepath

    filepath = '/home/kawaguchi/' // trim(dir_name)

    ! 配列初期化
    allocate(typ(1)%mol(NUM_MOL(1)))
    allocate(typ(2)%mol(NUM_MOL(2)))
    allocate(typ(3)%mol(NUM_MOL(3)))

        !open(6,*) は使用できない
        open(1,file=trim(filepath) // '/condition.dat')
    ! 各分子の位置データの出力
        open(10,file=trim(filepath) // '/posit_PtUp.dat')
        open(11,file=trim(filepath) // '/posit_Ar.dat')
        open(12,file=trim(filepath) // '/posit_PtLw.dat')
    ! 可視化用のpvch.fを移植 
        ! open(15,file=trim(filepath) // '/pos.dat')
        open(16,file=trim(filepath) // '/ovito.dat')
    ! 各分子の速度データの出力
        open(20,file=trim(filepath) // '/veloc_PtUp.dat')
        open(21,file=trim(filepath) // '/veloc_Ar.dat')
        open(22,file=trim(filepath) // '/veloc_PtLw.dat')
    ! 系のエネルギーデータの出力
        open(30,file=trim(filepath) // '/energy_PtUp.dat')
        open(31,file=trim(filepath) // '/energy_Ar.dat')
        open(32,file=trim(filepath) // '/energy_PtLw.dat')
        open(35,file=trim(filepath) // '/energy_all.dat')
    ! 系の温度データの出力
        open(40,file=trim(filepath) // '/temp.dat')
        open(45,file=trim(filepath) // '/temp_Layer.dat')

    ! 熱流束のデータ
        open(60,file=trim(filepath) // '/heatflux.dat')
        open(61,file=trim(filepath) // '/pressure.dat')
        
        open(70,file=trim(filepath) // '/force_phantom.dat')
        open(71,file=trim(filepath) // '/force_interface.dat')

    ! 各分子の最終位置データの出力
        open(80,file=trim(filepath) // '/finpos.dat')

    !--------------------PVWIN用
    !　分子の色
        ! open(90,file=trim(filepath) // '/mask.dat')

    ! write(15,'(3I7)') moltype, NUM_MOL_ALL, ndat
    ! do i = 1,ndat
    !     do j = 1, int(NUM_MOL(1)/NUM_Z(1))
    !         write(90,'(I7)') 15      ! 白色
    !     end do
    !     do j = int(NUM_MOL(1)/NUM_Z(1)) + 1, NUM_MOL(1)
    !         write(90,'(I7)') 14      ! 赤色
    !     end do
    !     do j = 1, NUM_MOL(2)
    !         write(90,'(I7)') 7       ! 黄色
    !     end do
    !     do j = 1, int(NUM_MOL(3)/NUM_Z(3))
    !         write(90,'(I7)') 15      ! 白色
    !     end do
    !     do j = int(NUM_MOL(3)/NUM_Z(3)) + 1, NUM_MOL(3)
    !         write(90,'(I7)') 0       ! 青色
    !     end do
    ! end do
    !---------------------------
    
    ! condition.datに設計条件を記録
    write(1,*) ''
    write(1, '(A18, F7.4, A3)') 'Scaling Time : ', time_scaling, ' ns '
    write(1, '(A18, F7.4, A3)') 'Relaxation Time : ', time_relax, ' ns '
    write(1, '(A18, F7.4, A3)') 'Measure Time : ', time_measure, ' ns '
    write(1,*) ''
    write(1, '(A8, I5)') 'Ar : ', NUM_MOL(2)
    write(1, '(A8, I5, A2)') 'Pt : ', NUM_MOL(1), '*2'
    write(1, '(A8, I5)') 'total : ', NUM_MOL_ALL
    write(1,*) ''
    write(1, '(A23, F4.2)') 'Strength Coefficient : ', INTER_STG
    write(1, '(A23, F5.1)') 'Temperature Ar : ', TEMP_Ar
    write(1, '(A23, F5.1)') 'Temperature Pt Upper : ', TEMP_Langevin(1)
    write(1, '(A23, F5.1)') 'Temperature Pt Lower : ', TEMP_Langevin(3)
    write(1,*) ''

    ! ターミナルに表示
    write(6,*) ''
    write(6, '(A18, I8, A4)') 'Scaling Step : ', step_scaling, 'step'
    write(6, '(A18, I8, A4)') 'Relaxation Step : ', step_relax, 'step'
    write(6, '(A18, I8, A4)') 'Measure Step : ', step_measure, 'step'
    write(6,*) ''
    write(6,*) '----------------------- Scaling Step -----------------------'
    write(6, '(I8, A4)') step_scaling, 'step'
    write(6,*) ''
    
    step_now = 0
    temp_layer_Ar(:) = 0.0d0
    temp_layer_Pt(:,:) = 0.0d0

    call seting ! 各分子の初期位置，初期速度などの設定
    call ovito

    do i = 1, step_all
        step_now = i
        ! ターミナルに表示
        if(step_now == step_scaling+1) then
            write(6,*) ''
            write(6,*) '---------------------- Relaxation Step ----------------------'
            write(6,'(I8, A4)') step_relax, 'step'
            write(6,*) ''
        end if
        if(step_now == step_scaling + step_relax+1) then
            write(6,*) ''
            write(6,*) '----------------------- Measure Step ------------------------'
            write(6,'(I8, A4)') step_measure, 'step'
            write(6,*) ''

            heat_phantom(:) = 0.0d0 ! パラメータモジュールで初期化するとうまくいかなかったのでここで初期化
            heat_interface(:) = 0.0d0
        end if

        call cpu_time(time)

        ! ステップ数が500の倍数のとき
        if (mod(step_now,500) == 0) then
            write(6,'(3X, I9, 1X, A3, I10, A1, F15.7, A3)') step_now, ' / ', step_all, '(', time, ' s)'
        endif

        ! スケーリング
        if (step_now <= step_scaling .and. mod(step_now,100) == 0) then
            call scaling ! 系内の全分子の温度の補正
        endif

        call calc_force_potential ! 各分子に働く力，速度，位置
        if(step_now > step_scaling + step_relax) then
            call calc_heatflux ! 熱流束の計算
        end if
        call calc_integration ! 位置、速度、加速度の更新

        call bound ! 境界条件の付与

        ! ステップ数が100の倍数+1のとき
        if(mod(step_now, 100) == 1) then
            call record_pos_vel ! 位置，速度を記録
            call record_energy_temp ! エネルギー，温度を記録
            
            if(step_now <= 100000) then ! 最初から最後まで記録するとデータ量が膨大になるので、適当に途中で止める
                call ovito ! OVITOという可視化ソフト用に記録を取る
            end if
        end if

        if(step_now > step_scaling + step_relax) then
            call record_pressure_heatflux ! 圧力、熱流束を記録
        end if
    end do

    ! step_now == 1 (mod 100) で記録してきたが、最後は step_now == 0 でも記録
    call calc_force_potential
    call calc_heatflux
    call calc_integration
    call record_energy_temp
    call record_pressure_heatflux
    call record_final ! 最終状態を記録

    write(6,*) ''
    write(6,*) '----------------------- Completed... ------------------------'
    
    contains

    ! 各分子の初期位置，初期速度などの設定
    subroutine seting
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: num, i, j, k
        double precision :: coord(3) = 0.0d0 ! xyz座標 (coordinate)
        double precision :: offset(3)
        double precision :: ran, alpha, beta, cr
        double precision :: v(3)

        do i = 1, 3
            coef_force(i) = 24.0d0*EPS(i)/SIG(i)  ! 無次元なことに注意 
        end do
            coef_force(4) = INTER_STG*24.0d0*EPS(4)/SIG(4)
        
        bnd_len(:) = BND_LEN0(:)

        ! PVWIN用
        ! write(15,*)BND_LEN0(1), BND_LEN0(2), BND_LEN0(3)
        ! write(15,*)ntime0, ndt

        do i = 1, TYPMOL
            cutoff(:) = bnd_len(:) - CUTOFF*SIG(i)
        end do

        num = 0

        !上段のPt配置
        offset(1) = ST_DIST(1)*0.25d0
        offset(2) = ST_DIST(1)*0.25d0
        offset(3) = BND_LEN_Z0 - ST_DIST(1)*0.25d0
        do k = 1,NUM_Z(1)
            coord(3) = offset(3) - dble(k-1)*ST_DIST(1)*0.5d0
            do i = 1,NUM_X(1)
                coord(1) = offset(1) + dble(i-1)*ST_DIST(1)*0.5d0
                do j = 1,NUM_Y(1)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(1)   !x偶数
                        else
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(1) + ST_DIST(1)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(1) + ST_DIST(1)*0.5d0    !x偶数
                        else
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(1)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(1)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        !中段のAr配置
        num = 0
        offset(1) = BND_LEN0(1)*0.50d0 - ST_DIST(2)*(0.25d0*(NUM_X(2)  -1))
        offset(2) = BND_LEN0(2)*0.50d0 - ST_DIST(2)*(0.25d0*(NUM_Y(2)*2-1))
        offset(3) = BND_LEN0(3)*0.50d0 - ST_DIST(2)*(0.25d0*(NUM_Z(2)  -1))
        do k = 1,NUM_Z(2)
            coord(3) = offset(3) + dble(k-1)*ST_DIST(2)*0.5d0
            do i = 1,NUM_X(2)
                coord(1) = offset(1) + dble(i-1)*ST_DIST(2)*0.5d0
                do j = 1,NUM_Y(2)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(2)   !x偶数
                        else
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(2) + ST_DIST(2)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(2) + ST_DIST(2)*0.5d0    !x偶数
                        else
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(2)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(2)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        !下段のPt配置
        num = 0
        offset(1) = ST_DIST(3)*0.25d0
        offset(2) = ST_DIST(3)*0.25d0
        offset(3) = ST_DIST(3)*0.25d0
        do k = 1,NUM_Z(3)
            coord(3) = offset(3) + dble(k-1)*ST_DIST(3)*0.5d0
            do i = 1,NUM_X(3)
                coord(1) = offset(1) + dble(i-1)*ST_DIST(3)*0.5d0
                do j = 1,NUM_Y(3)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(3)   !x偶数
                        else
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(3) + ST_DIST(3)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(3) + ST_DIST(3)*0.5d0    !x偶数
                        else
                            coord(2) = offset(2) + dble(j-1)*ST_DIST(3)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(3)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        cr = 1.0d-6
        do j = 1, TYPMOL
            do i = 1, NUM_MOL(j)
                call random_number(ran)
                alpha = PI*ran
                call random_number(ran)
                beta = 2.000d0*PI*ran
                v(1) = dsin(alpha)*dcos(beta)*cr
                v(2) = dsin(alpha)*dsin(beta)*cr
                v(3) = dcos(alpha)*cr
                typ(j)%mol(i)%vel(:) = v(:)
            end do
        end do

        do j = 1, TYPMOL
            if(j == 2) then
                cycle
            end if

            do i = 1, NUM_X(j)*NUM_Y(j)        
                typ(j)%mol(i)%vel(:) = 0.0d0
            end do
        end do
    end subroutine seting

    ! 系内の全分子の温度の補正
    subroutine scaling
        use parameters
        use variable
        use molecules_struct
        implicit none
        double precision :: temptp, vel2, aimtem, aimnot, baiss
        integer :: i

        temptp = 0.0d0
        do i = 1, NUM_MOL(2)
            vel2 = typ(2)%mol(i)%vel(1)**2 + typ(2)%mol(i)%vel(2)**2 + typ(2)%mol(i)%vel(3)**2
            temptp = temptp + vel2
        end do
        temptp = temptp / NUM_MOL(2) * 1.0d-16
        aimtem = TEMP_Ar
        aimnot = 3.000d0 * BOLTZ * aimtem / MASS(2)
        baiss = dsqrt(aimnot / temptp)

        ! 速度ベクトルのスケーリング
        ! Arのみ
        do i = 1, NUM_MOL(2)
            typ(2)%mol(i)%vel(:) = typ(2)%mol(i)%vel(:) * baiss
        end do
    end subroutine scaling

    subroutine ovito
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        double precision :: t
        
        t = dble(step_now)*dt*1.0d-6 * 100 ! 100ステップおきに計測　単位は [ns]

        write(16, *) NUM_MOL_ALL
        write(16, '(A8, F4.1, A13, F4.1, A13, F4.1, A39, F6.3)') 'Lattice=', BND_LEN_X0, '0.0 0.0 0.0 ', &
                BND_LEN_Y0, '0.0 0.0 0.0 ', BND_LEN_Z0, '"Properties="species:S1:pos:R:3" Time=', t
        do j = 1, TYPMOL
            if(j == 2) then
                do i = 1, NUM_MOL(j)
                    write(16, '(A2, 3E15.7)') 'Ar', typ(j)%mol(i)%pos(1), typ(j)%mol(i)%pos(2), typ(j)%mol(i)%pos(3)
                end do
            else
                do i = 1, NUM_MOL(j)
                    write(16, '(A2, 3E15.7)') 'Pt', typ(j)%mol(i)%pos(1), typ(j)%mol(i)%pos(2), typ(j)%mol(i)%pos(3)
                end do
            end if
        end do
    end subroutine ovito

    subroutine calc_force_potential
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j, k, i1, i2

        integer, parameter :: SAFE = 2 ! 安全率
        double precision :: rbook(5000, TYPMOL) ! 5000は適当に決めました笑
        logical :: judge_same(5000, 5000, TYPMOL), judge_diff(NUM_MOL(1), NUM_MOL(2), TYPMOL) ! 相互作用を計算する分子を区別するためのフラグ
        
        double precision :: rnd
        double precision :: div(3), dist
        double precision :: dit2, dit4, dit6, dit8, dit12, dit14
        double precision :: ppp, force, forVec(3)
        double precision :: vel, vmax(5000, TYPMOL), vene(3), sumvene
    
        do j = 1, TYPMOL
            do i = 1, NUM_MOL(j)
                typ(j)%mol(i)%acc(:) = 0.0d0
                typ(j)%mol(i)%pot  = 0.0d0
                typ(j)%mol(i)%kin  = 0.0d0
            end do
        end do
        force_inter(:,:,:) = 0.0d0
    
        ! 粒子登録法で使うための速度
        do j = 1, TYPMOL
            do i = 1, NUM_MOL(j)
                vel = dsqrt(typ(j)%mol(i)%vel(1)**2 + typ(j)%mol(i)%vel(2)**2 + typ(j)%mol(i)%vel(3)**2)
                if(vmax(i,j) < vel) then
                    vmax(i,j) = vel
                end if
            end do
        end do

        ! 同じ分子同士の影響
            ! 粒子登録法
            if(mod(step_now,step_update) == 1) then ! 最大速度を更新
                do j = 1, TYPMOL
                    do i = 1, NUM_MOL(j)
                        rbook(i,j) = CUTOFF*SIG(j) + SAFE*step_update*vmax(i,j)*dt
                        vmax(i,j) = 0.0d0
                    end do

                    do i1 = 1, NUM_MOL(j)
                        do i2 = i1+1, NUM_MOL(j)
                            div(:) = typ(j)%mol(i1)%pos(:) - typ(j)%mol(i2)%pos(:)
    
                            do k = 1, 2
                                if(div(k) < -bnd_len(k)/2) then
                                    div(k) = div(k) + bnd_len(k)
                                else if(div(k) > bnd_len(k)/2) then
                                    div(k) = div(k) - bnd_len(k)
                                end if
                            end do

                            dit2 = div(1)**2 + div(2)**2 + div(3)**2
                            dist = dsqrt(dit2)

                            if(rbook(i1,j) > dist .or. rbook(i2,j) > dist) then
                                judge_same(i1,i2,j) = .true.
                            else
                                judge_same(i1,i2,j) = .false.
                            end if
                        end do
                    end do
                end do
            end if

            ! 分子間の相互作用力 → ポテンシャルエネルギー
            do j = 1, TYPMOL
                do i1 = 1, NUM_MOL(j)
                    do i2 = i1+1, NUM_MOL(j)
                        if(.not. judge_same(i1,i2,j)) then
                            cycle
                        else
                            div(:) = typ(j)%mol(i1)%pos(:) - typ(j)%mol(i2)%pos(:)
                            ! カットオフ
                            do k = 1, 3
                                if (div(k) < -cutoff(k)) then
                                    div(k) = div(k) + bnd_len(k)
                                else if(div(k) > cutoff(k)) then
                                    div(k) = div(k) - bnd_len(k)
                                endif
                
                                div(k) = div(k) / SIG(j)
            
                                if (abs(div(k)) > CUTOFF) then
                                    cycle
                                endif
                            end do
            
                            dit2 = div(1)**2 + div(2)**2 + div(3)**2
                            dist = dsqrt(dit2)

                            if(dist > CUTOFF) then
                                cycle
                            endif
            
                            dit4   = dit2*dit2
                            dit6   = dit4*dit2
                            dit8   = dit4*dit4
                            dit12  = dit6*dit6
                            dit14  = dit8*dit6
                            ppp    = 4.00d0*EPS(j)*(1.0d0/dit12-1.0d0/dit6)
                            typ(j)%mol(i1)%pot = typ(j)%mol(i1)%pot + ppp*0.500d0
                            typ(j)%mol(i2)%pot = typ(j)%mol(i2)%pot + ppp*0.500d0
            
                            force  = coef_force(j)*(-2.00d0/dit14+1.0d0/dit8)
                            forVec(:) = -force*div(:)
                            typ(j)%mol(i1)%acc(:) = typ(j)%mol(i1)%acc(:) + forVec(:)/MASS(j)
                            typ(j)%mol(i2)%acc(:) = typ(j)%mol(i2)%acc(:) - forVec(:)/MASS(j)
                        end if
                    end do
                end do
            end do
    
        ! 異なる分子同士の影響
            ! 粒子登録法
            if(mod(step_now,step_update) == 1) then
                do j = 1, TYPMOL
                    if(j == 2) then
                        cycle
                    end if
                    do i1 = 1, NUM_MOL(j) ! Pt
                        do i2 = 1, NUM_MOL(2) ! Ar
                            div(:) = typ(j)%mol(i1)%pos(:) - typ(2)%mol(i2)%pos(:)
    
                            do k = 1, 2
                                if(div(k) < -bnd_len(k)/2) then
                                    div(k) = div(k) + bnd_len(k)
                                else if(div(k) > bnd_len(k)/2) then
                                    div(k) = div(k) - bnd_len(k)
                                end if
                            end do

                            dit2 = div(1)**2 + div(2)**2 + div(3)**2
                            dist = dsqrt(dit2)

                            if(rbook(i2,2) > dist) then     ! Arの条件のみでよい？
                                judge_diff(i1,i2,j) = .true.
                            else
                                judge_diff(i1,i2,j) = .false.
                            end if
                        end do
                    end do
                end do
            end if

            ! Ar-Ptの処理    配列は1が上Pt, 2がAr, 3がしたPt
            do j = 1, TYPMOL
                if(j == 2) then
                    cycle
                end if
                do i1 = 1, NUM_MOL(j)       ! Pt
                    do i2 = 1, NUM_MOL(2)    ! Ar
                        if(.not. judge_diff(i1,i2,j)) then
                            cycle
                        else
                            div(:) = typ(j)%mol(i1)%pos(:) - typ(2)%mol(i2)%pos(:)
                
                            do k = 1, 3
                                if (div(k) < -cutoff(k)) then
                                    div(k) = div(k) + bnd_len(k)
                                else if(div(k) > cutoff(k)) then
                                    div(k) = div(k) - bnd_len(k)
                                endif
                
                                div(k) = div(k) / SIG(4)
                
                                if (abs(div(k)) > CUTOFF) then
                                    cycle
                                endif
                            end do
                
                            dit2 = div(1)**2 + div(2)**2 + div(3)**2
                            dist = dsqrt(dit2)
                
                            if(dist > CUTOFF) then
                                cycle
                            endif
                
                            dit4   = dit2*dit2
                            dit6   = dit4*dit2
                            dit8   = dit4*dit4
                            dit12  = dit6*dit6
                            dit14  = dit8*dit6
                            ppp    = INTER_STG*4.00d0*EPS(4)*(1.0d0/dit12-1.0d0/dit6) ! 異分子間ではINTER_STG(接触角)を忘れずに
                            typ(j)%mol(i1)%pot = typ(j)%mol(i1)%pot + ppp*0.500d0
                            typ(2)%mol(i2)%pot = typ(2)%mol(i2)%pot + ppp*0.500d0
                
                            force  = coef_force(4)*(-2.00d0/dit14+1.0d0/dit8)
                            forVec(:) = -force*div(:)
                            force_inter(i1,j,:) = force_inter(i1,j,:) - forVec(:) ! 無次元なことに注意　符号が逆な気がする
                
                            typ(j)%mol(i1)%acc(:) = typ(j)%mol(i1)%acc(:) + forVec(:)/MASS(j)
                            typ(2)%mol(i2)%acc(:) = typ(2)%mol(i2)%acc(:) - forVec(:)/MASS(2)
                        end if
                    end do
                end do
            end do
        
        ! PtのPhantom層はダンパー力とランダム力を付与
        do j = 1, TYPMOL
            if(j == 2) then
                cycle
            end if
    
            do i = NUM_X(j)*NUM_Y(j) + 1, 2*NUM_X(j)*NUM_Y(j) ! Phantom層のみ
                do k = 1, 3
                    rnd = Random() 
                    ! ランダム力
                    force_rand(i,j,k) = rnd * getStddev(TEMP_Langevin(j)) ! 標準偏差を無次元化したまま扱う
                    ! ダンパー力
                    force_damp(i,j,k) = -DAMP * typ(j)%mol(i)%vel(k) * 1.0d+14 ! ダンパー力を無次元化したまま扱う
                end do
                ! ランダム力とダンパー力を追加
                typ(j)%mol(i)%acc(:) = typ(j)%mol(i)%acc(:) + (force_rand(i,j,:) + force_damp(i,j,:)) / MASS(j)*1.0d-3 ! -9+26-20 = -3  (ランダム, ダンパー力の有次元化-9)*(質量の有次元化+26)*(加速度の無次元化-20)
            end do
        end do

        ! 運動エネルギー計算
        do j = 1, TYPMOL
            do i = 1, NUM_MOL(j)
                typ(j)%mol(i)%vtmp(:) = typ(j)%mol(i)%vel(:) + typ(j)%mol(i)%acc(:)*0.500d0*dt ! vel(t) = vel(t-dt/2) + acc(t)*dt/2
                sumvene = typ(j)%mol(i)%vtmp(1)**2 + typ(j)%mol(i)%vtmp(2)**2 + typ(j)%mol(i)%vtmp(3)**2
                typ(j)%mol(i)%kin = 0.500d0*MASS(j)*sumvene
            end do
        end do   
    end subroutine calc_force_potential

    subroutine calc_heatflux
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j

        do j = 1, TYPMOL
            if (j == 2) then
                cycle
            end if
    
            do i = NUM_X(j)*NUM_Y(j) + 1, 2*NUM_X(j)*NUM_Y(j) ! Phantom層  
                heat_phantom(j) = heat_phantom(j) + &
                    &  (force_rand(i,j,1) + force_damp(i,j,1)) * typ(j)%mol(i)%vtmp(1) + &
                    &  (force_rand(i,j,2) + force_damp(i,j,2)) * typ(j)%mol(i)%vtmp(2) + &
                    &  (force_rand(i,j,3) + force_damp(i,j,3)) * typ(j)%mol(i)%vtmp(3) ! 速さの有次元化 10^5
            end do
            
            do i = 1, NUM_MOL(j) ! Pt分子全体
                heat_interface(j) = heat_interface(j) + &
                    &  force_inter(i,j,1) * typ(j)%mol(i)%vtmp(1) + &
                    &  force_inter(i,j,2) * typ(j)%mol(i)%vtmp(2) + &
                    &  force_inter(i,j,3) * typ(j)%mol(i)%vtmp(3)
            end do
        end do
    end subroutine calc_heatflux

    subroutine calc_integration
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j

        ! 数値積分 (蛙跳び法)
        do j = 1, TYPMOL
            ! Arの計算
            if(j == 2) then
                do i = 1, NUM_MOL(2)
                    typ(2)%mol(i)%vel(:) = typ(2)%mol(i)%vel(:) + typ(2)%mol(i)%acc(:) * dt   ! vel(t+dt/2) = vel(t-dt/2) + acc(t)*dt
                    typ(2)%mol(i)%pos(:) = typ(2)%mol(i)%pos(:) + typ(2)%mol(i)%vel(:) * dt   ! pos(t+dt)   = pos(t)      + vel(t+dt/2)*dt
                end do
            else
            ! Ptの計算
                ! 固定層
                do i = 1, NUM_X(j)*NUM_Y(j)
                    typ(j)%mol(i)%vel(:) = 0.0d0
                end do
    
                ! その他の層
                do i = NUM_X(j)*NUM_Y(j) + 1, int(NUM_MOL(j))
                    typ(j)%mol(i)%vel(:) = typ(j)%mol(i)%vel(:) + typ(j)%mol(i)%acc(:) * dt
                    typ(j)%mol(i)%pos(:) = typ(j)%mol(i)%pos(:) + typ(j)%mol(i)%vel(:) * dt
                end do
            end if
        end do
    end subroutine calc_integration

    subroutine bound ! 境界条件の付与
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        do j = 1, TYPMOL
            do i = 1, NUM_MOL(j)
                if(typ(j)%mol(i)%pos(1) < 0.00d0) then
                    typ(j)%mol(i)%pos(1) = typ(j)%mol(i)%pos(1) + bnd_len(1)
                else if(typ(j)%mol(i)%pos(1) > bnd_len(1)) then
                    typ(j)%mol(i)%pos(1) = typ(j)%mol(i)%pos(1) - bnd_len(1)
                endif

                if(typ(j)%mol(i)%pos(2) < 0.00d0) then
                    typ(j)%mol(i)%pos(2) = typ(j)%mol(i)%pos(2) + bnd_len(2)
                else if(typ(j)%mol(i)%pos(2) > bnd_len(2)) then
                    typ(j)%mol(i)%pos(2) = typ(j)%mol(i)%pos(2) - bnd_len(2)
                endif

                ! z方向は周期境界の補正を行わない
                ! if(typ(j)%mol(i)%pos(3) < 0.00d0) then
                !     typ(j)%mol(i)%pos(3) = typ(j)%mol(i)%pos(3) + bnd_len(3)
                ! else if(typ(j)%mol(i)%pos(3) > bnd_len(3)) then
                !     typ(j)%mol(i)%pos(3) = typ(j)%mol(i)%pos(3) - bnd_len(3)
                ! endif
            end do
        end do
    end subroutine bound

    subroutine record_pos_vel ! 位置，速度を記録
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j, k
    
        do i = 1, NUM_MOL(1)
            ! posit_PtUp.dat
            write(10, '(I6, 3E15.7)') (i, typ(1)%mol(i)%pos(k), k = 1,3)
            ! veloc_PtUp.dat
            write(20, '(I6, 3E15.7)') (i, typ(1)%mol(i)%vel(k), k = 1,3)
        end do
        do i = 1, NUM_MOL(2)
            ! posit_Ar.dat
            write(11, '(I6, 3E15.7)') (i, typ(2)%mol(i)%pos(k), k = 1,3)
            ! veloc_Ar.dat
            write(21, '(I6, 3E15.7)') (i, typ(2)%mol(i)%vel(k), k = 1,3)
        end do
        do i = 1, NUM_MOL(3)
            ! posit_PtLw.dat
            write(12, '(I6, 3E15.7)') (i, typ(3)%mol(i)%pos(k), k = 1,3)
            ! veloc_PtLw.dat
            write(22, '(I6, 3E15.7)') (i, typ(3)%mol(i)%vel(k), k = 1,3)
        end do
    
        do i = NUM_X(1)*NUM_Y(1) + 1, 2*NUM_X(1)*NUM_Y(1) ! Phantom層
            write(70, '(I6, 6E15.7)') i, force_rand(i,1,1)*1.0d-9, force_rand(i,1,2)*1.0d-9, force_rand(i,1,3)*1.0d-9, &
                                         force_damp(i,1,1)*1.0d-9, force_damp(i,1,2)*1.0d-9, force_damp(i,1,3)*1.0d-9
        end do
    
        do i = 1, NUM_MOL(1)!(NUM_Z(1)-1)*NUM_X(1)*NUM_Y(1) + 1, NUM_MOL(1) ! 固液界面層
            write(71, '(I6, 3E15.7)') i, force_inter(i,1,1)*1.0d-6, force_inter(i,1,2)*1.0d-6, force_inter(i,1,3)*1.0d-6 ! 上Pt, Ar, 下Pt
        end do
    
        ! PVWIN用
        ! do j = 1, TYPMOL
        !     do i = 1, NUM_MOL(j)
        !         ! pos.dat
        !         write(15, '(3E15.7)') typ(j)%mol(i)%pos(1), typ(j)%mol(i)%pos(2), typ(j)%mol(i)%pos(3)
        !     end do
        end do
    end subroutine record_pos_vel
    
    subroutine record_energy_temp ! エネルギー，温度を記録
        use parameters
        use variable
        use molecules_struct
        implicit none
        double precision, dimension(TYPMOL) :: ene_sum, pot_sum, kin_sum, temp
        double precision :: ene_all, pot_all, kin_all
        double precision :: kin_Pt_tmp, kin_Ar_tmp(NUM_DIV_Ar)
        integer :: cnt(NUM_DIV_Ar)
        double precision :: temp_layer_Pt_(NUM_Z(1),TYPMOL)
        double precision :: temp_layer_Ar_(NUM_DIV_Ar)
        integer :: i, j, k
        integer :: step
        step = step_now-step_relax-step_scaling
    
        ! sum は種類ごとの総和、all はArとPtを合わせた全分子の総和で一応分けている
        ene_all = 0.0d0
        pot_all = 0.0d0
        kin_all = 0.0d0
        ene_sum(:) = 0.0d0
        pot_sum(:) = 0.0d0
        kin_sum(:) = 0.0d0
        temp(:) = 0.0d0
        temp_layer_Pt_(:,:) = 0.0d0
        temp_layer_Ar_(:) = 0.0d0
    
        do j = 1, TYPMOL
            ! ポテンシャル
            do i = 1, NUM_MOL(j)
                pot_sum(j) = pot_sum(j) + typ(j)%mol(i)%pot
            end do
            pot_sum(j) = pot_sum(j) * 1.0d-16
    
            ! 運動エネルギー
            if (j == 2) then
            ! Arの温度
                kin_Ar_tmp(:) = 0.0d0 ! 各領域のAr分子の運動エネルギーを足すための変数
                cnt(:) = 0
                do i = 1, NUM_MOL(j)
                    do k = 1, NUM_DIV_Ar
                        if ( (k-1)*LEN_DIV_Ar <= (typ(j)%mol(i)%pos(3) - THICK(3)) .and. (typ(j)%mol(i)%pos(3) - THICK(3)) < k*LEN_DIV_Ar) then
                            kin_Ar_tmp(k) = kin_Ar_tmp(k) + typ(j)%mol(i)%kin
                            cnt(k) = cnt(k) + 1
                            cycle
                        end if
                    end do
                end do
    
                do k = 1, NUM_DIV_Ar
                    kin_Ar_tmp(k) = kin_Ar_tmp(k) * 1.0d-16
                    kin_sum(j) = kin_sum(j) + kin_Ar_tmp(k)
                    temp_layer_Ar_(k) = 2.000d0 * kin_Ar_tmp(k) / (3.000d0 * dble(cnt(k)) * BOLTZ)
                    temp(j) = temp(j) + temp_layer_Ar_(k)
                end do
    
                temp(j) = temp(j) / dble(NUM_DIV_Ar)
            else
            ! Ptの温度
                do k = 2, NUM_Z(j) ! Ptの層の数
                    kin_Pt_tmp = 0.0d0
                    do i = (k-1)*NUM_X(j)*NUM_Y(j) + 1, k*NUM_X(j)*NUM_Y(j)
                        kin_Pt_tmp = kin_Pt_tmp + typ(j)%mol(i)%kin
                    end do
    
                    kin_Pt_tmp = kin_Pt_tmp * 1.0d-16
                    kin_sum(j) = kin_sum(j) + kin_Pt_tmp
                    temp_layer_Pt_(k, j) = 2.000d0 * kin_Pt_tmp / (3.000d0 * dble(NUM_MOL(j)/NUM_Z(j)) * BOLTZ)
                    temp(j) = temp(j) + temp_layer_Pt_(k, j)
                end do
    
                temp(j) = temp(j) / dble(NUM_Z(j)-1)
            end if
    
            ene_sum(j) = pot_sum(j) + kin_sum(j)
            ene_all = ene_all + ene_sum(j)
            pot_all = pot_all + pot_sum(j)
            kin_all = kin_all + kin_sum(j)
        end do
    
        if(step_now > step_scaling + step_relax) then
            do j = 1, TYPMOL
                if(j == 2) then
                    do i = 1, NUM_DIV_Ar
                        temp_layer_Ar(i) = temp_layer_Ar(i) + temp_layer_Ar_(i)
                    end do
                else
                    do i = 2, NUM_Z(j)
                        temp_layer_Pt(i,j) = temp_layer_Pt(i,j) + temp_layer_Pt_(i,j)
                    end do
                end if
            end do
        end if
        
    
        write(30, '(5E15.7)') step_now*int(dt)*1.0d-6, ene_sum(1), pot_sum(1), kin_sum(1)  ! energy_PtUp.dat
        write(31, '(5E15.7)') step_now*int(dt)*1.0d-6, ene_sum(2), pot_sum(2), kin_sum(2)  ! energy_Ar.dat
        write(32, '(5E15.7)') step_now*int(dt)*1.0d-6, ene_sum(3), pot_sum(3), kin_sum(3)  ! energy_PtLw.dat
        write(35, '(5E15.7)') step_now*int(dt)*1.0d-6, ene_all, pot_all, kin_all           ! energy_all.dat
        write(40, '(5E15.7)') step_now*int(dt)*1.0d-6, temp(1), temp(2), temp(3)        ! temp.dat
    end subroutine record_energy_temp
    
    subroutine record_pressure_heatflux ! 熱流束を記録
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        double precision :: heat_phantom_sum(TYPMOL) ! Phantom層からの熱輸送量
        double precision :: heat_interface_sum(TYPMOL) ! 固液界面での熱輸送量
        double precision :: pressure(TYPMOL) ! 圧力
        integer :: step
        step = step_now-step_relax-step_scaling
    
        heat_phantom_sum(:) = 0.0d0
        heat_interface_sum(:) = 0.0d0
        pressure(:) = 0.0d0
    
        do j = 1, TYPMOL
            if (j == 2) then
                cycle
            end if
    
            heat_phantom_sum(j) = heat_phantom(j) / AREA_Pt * 1.0d+1 ! -9+5+20-15 = +1  (ランダム力-9)*(速度+5)*(面積の逆数+20)*(時間-15)
            heat_interface_sum(j) = heat_interface(j) / AREA_Pt * 1.0d+4 ! -6+5+20-15 = +4  (相互作用力-6)*(速度+5)*(面積の逆数+20)*(時間-15)
            
            do i = 1, NUM_MOL(j)
                if(j == 1) then
                    pressure(1) = pressure(1) + force_inter(i,1,3) / AREA_Pt * 1.0d+8 ! [MPa]　-6+20-6 = +8  (相互作用力-6)*(面積の逆数+20)*(メガ-6)
                else
                    pressure(3) = pressure(3) - force_inter(i,3,3) / AREA_Pt * 1.0d+8
                end if
            end do
        end do
    
        write(60,'(5E15.7)') step*int(dt)*1.0d-6, heat_phantom_sum(1), heat_phantom_sum(3), heat_interface_sum(1), heat_interface_sum(3)
        write(61,'(4E15.7)') step*int(dt)*1.0d-6, pressure(1), pressure(2), pressure(3)
    end subroutine record_pressure_heatflux
    
    subroutine record_final ! 最終状態を記録
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        double precision :: disz
    
        do j = 1, TYPMOL
            if(j == 2) then
                do i = 1, NUM_DIV_Ar
                    temp_layer_Ar(i) = temp_layer_Ar(i) / dble(step_measure/100)
                end do
            else
                do i = 2, NUM_Z(j)
                    temp_layer_Pt(i,j) = temp_layer_Pt(i,j) / dble(step_measure/100)
                end do
            end if
        end do
    
        disz = ST_DIST(3)*0.25d0
        do i = 2, NUM_Z(3)   ! 固定層は除外
            disz = disz + ST_DIST(3)*0.5d0
            write(45, '(2E15.7)') disz, temp_layer_Pt(i,3)
        end do
    
        disz = disz - LEN_DIV_Ar * 0.5d0
        do i = 1, NUM_DIV_Ar
            disz = disz + LEN_DIV_Ar
            write(45, '(2E15.7)') disz, temp_layer_Ar(i)
        end do
    
        disz = BND_LEN_Z0 - ST_DIST(1)*0.25d0
        do i = 2, NUM_Z(1)
            disz = disz - ST_DIST(1)*0.5d0
            write(45, '(2E15.7)') disz, temp_layer_Pt(i,1)
        end do
    
    end subroutine record_final
end program main

!--------追記--------!
! 繰り返し処理で j == 2のとき cycle（スキップ）という処理が多いと思いますが、これは分子の種類の順番を上Pt, Ar, 下Ptとしてる影響です
! いちいちめんどくさいので、上Pt, 下Pt, Arという順番にしましょう

! 変数名の重複を避けるために、あまり使わない変数は一つのモジュールで完結するよう心がけました
! メインプログラムの call の一行のみをコメントアウトすることで簡単に処理のオンオフができるようにしました（例えば位置、速度の記録は動作の確認でしか使わないので、後半はずっとコメントアウトしていました）
! モジュールの分け方は、たくさん分けすぎるとオフしにくくなるので、一括りの処理でまとめるのをおすすめします（最初はcalc~ をすべてcalcというモジュールにまとめていましたが、あまりに長いので分けました）

! 複数のファイルにコードを分けるのもあまりおすすめしません
! 置換が面倒だったり、実行コマンドが長くなります
! もし分けるとしても、変数定義ファイルとプログラムファイルの二つに分けるのが一般的だと思います

! 最初の方に構造体を作ってみて計算時間にあまり影響がなかったので未だに使い続けていますが、コードが長くなるので pos(kind, number) という感じで変数定義したほうがいいです
! さらに多次元配列にすると計算時間が長くなる可能性があるのでやめましょう