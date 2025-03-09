module parameters
    implicit none
!--------------- よく変更する---------------------!
    ! 計算時間
    double precision, parameter :: timeScaling = 0.1d0 ! スケーリング時間[ns]
    double precision, parameter :: timeRelax   = 1.0d0 ! 緩和計算時間[ns]
    double precision, parameter :: timeMeasure = 10.0d0 ! データ計測時間[ns]
    double precision, parameter :: timeSum = timeScaling + timeRelax + timeMeasure ! 全計算時間
    
    double precision, parameter :: dt = 1.00d0 ! 無次元時間ステップ(無次元，有次元の場合fs)
    double precision, parameter :: tau = 1.00d0 ! 測定間隔[fs]

    double precision, parameter :: tempAr = 100d0 ! 系内（目標）温度  K
    double precision, parameter :: angCon = 0.01d0 ! 相互作用強さ

    character(len=20) :: dir_name = 'a0.01_'  ! ファイル名

    ! ステップ
    integer, parameter :: stpScaling = int(timeScaling/dt*1.0d+6) ! スケーリングステップ
    integer, parameter :: stpRelax   = int(timeRelax  /dt*1.0d+6) ! 緩和計算ステップ
    integer, parameter :: stpMeasure = int(timeMeasure/dt*1.0d+6) ! データ計測ステップ
    integer, parameter :: stpMax = stpScaling + stpRelax + stpMeasure ! 最大ステップ数

    ! 分子の数，配置
    integer, parameter :: TYPMOL = 3 ! 分子の種類数
    integer, parameter :: COMP = 3 ! 相互作用の組み合わせ数 = nC2
    integer, parameter :: numx(TYPMOL) = [20,10,20]
    integer, parameter :: numy(TYPMOL) = [10, 8,10]
    integer, parameter :: numz(TYPMOL) = [ 6,25, 6]
    integer, parameter :: nummol(TYPMOL) = [numx(1)*numy(1)*numz(1), numx(2)*numy(2)*numz(2), numx(3)*numy(3)*numz(3)] ! 各分子の数
    integer, parameter :: tlnkoss = nummol(1) + nummol(2) + nummol(3)
    
    ! 境界条件
    double precision, parameter :: STDIST(TYPMOL) = [3.92d0, 4.8d0, 3.92d0] ! 格子定数(無次元)[Å]
    double precision, parameter :: thick1 = STDIST(1)*(numz(1)*0.5d0+0.25d0)
    double precision, parameter :: thick2 = STDIST(2)*numz(2)*0.5d0
    double precision, parameter :: thick3 = STDIST(3)*(numz(3)*0.5d0+0.25d0)
    double precision, parameter :: thick(TYPMOL) = [thick1, thick2, thick3]
    double precision, parameter :: xsyul0 = STDIST(1) * numy(1) ! x方向の周期境界長さ(無次元)
    double precision, parameter :: ysyul0 = STDIST(1) * numy(1) ! y方向の周期境界長さ(無次元）
    double precision, parameter :: zsyul0 = 87.0d0 ! thick(1)+thick(2)+thick(3) ! z方向の周期境界長さ(無次元）
    double precision, parameter :: syul0(3) = [xsyul0, ysyul0, zsyul0]
    
    double precision, parameter :: CUTOFF = 3.300d0 ! カットオフ長さ/σ
    double precision, parameter :: AVOGA = 6.022d+23 ! アボガドロ数
    double precision, parameter :: BOLTZ = 1.3806662d-23 ! ボルツマン定数 [J/K]
    double precision, parameter :: PI = 3.141592654d0 ! 円周率

    ! 分子の質量
    !double precision, parameter :: bunsi(TYPMOL) = [195.084d-3, 39.950d-3, 195.084d-3] ! 分子の質量  kg/mol   
    double precision, parameter :: MASS(TYPMOL) = [32.395d0, 6.6340d0, 32.395d0] ! 分子の質量（無次元） * 10d-26 [kg/個]
    ! Lennard-Jonesパラメータ
    double precision, parameter :: SIG(COMP+1) = [2.540d0, 3.400d0, 2.540d0, 2.970d0]  ! σ(無次元) *1.0d-10
    double precision, parameter :: EPS(COMP+1) = [109.2d-5, 1.666d-5, 109.2d-5, 13.49d-5] ! ε(無次元) *1.0d-16

    ! 粒子登録法
    integer, parameter :: stpUpdate = 40 ! 更新ステップ

    ! Langevin法
    double precision, parameter :: tempLanPt(TYPMOL) = [130d0, 0d0, 70d0] ! Langevin法を用いるPtの温度  真ん中は使わない
    double precision, parameter :: DIRAC = 1.054571817d-34 ! ディラック定数 [J･s]
    double precision, parameter :: DEBTMP = 240d0 ! Debye温度 [K]
    double precision, parameter :: OMEGA = 3.14212728482d+13 ! BOLTZ * DEBTMP / DIRAC * 1.0d-11 ! Debye定数 (有次元)
    double precision, parameter :: DAMP =  5.32967075080d-12 ! MASS(1) * PI * OMEGA / 6.000d0 ! ダンパーの減衰係数 (有次元)
    ! 熱流束
    double precision, parameter :: areaPt = STDIST(1)*STDIST(1)*int(numx(1)*0.5)*numy(1)
    
    integer, parameter :: numDivAr = 15 ! Arの温度分布の分割数

    integer, parameter :: stpRecord = 1 ! 記録するステップ数間隔

end module parameters

! 変数
module variable
    use parameters
    implicit none
    integer :: stpNow ! 現在のステップ数

    ! 分子間力
    double precision :: forCoef(4) = [0.0d0, 0.0d0, 0.0d0, 0.0d0] ! 力の計算の係数

    ! カットオフ
    double precision :: cutof(3) ! ポテンシャルのカットオフ長さx,y,z方向
    double precision :: syul(3) ! ポテンシャルのカットオフ長さx,y,z方向，x,y,z方向の周期長さ

    ! Langevin法
    double precision :: rndForce(nummol(1),TYPMOL,3)   ! ランダム力用
    double precision :: dmpForce(nummol(1),TYPMOL,3)   ! ダンパー力用
    double precision :: stddev ! 標準偏差
    double precision :: rnd_2
    logical :: isOdd = .true. ! 乱数のsinとcosを交互に出すためのフラグ

    ! 熱流束
    double precision :: interForce(nummol(1),TYPMOL,3) ! 熱流束を計算するための相互作用力 z方向のみを使う
    double precision :: heatPhantom(TYPMOL)! = 0.0d0 ! Phantom層からの熱輸送量
    double precision :: heatInterface(TYPMOL)! = 0.0d0 ! 固液界面での熱輸送量

    double precision :: tempLayerPt(numz(1),TYPMOL) ! Ptの層ごとの温度
    double precision :: tempLayerAr(numDivAr) ! z方向に分割した領域内のArの温度
    double precision :: zdiv = thick(2)/numDivAr ! Arの温度分布の分割距離

    contains

    function Random() result(rnd_) ! 乱数生成用の関数　　呼び出す時に何回も計算しなくていいように下の関数と分けた
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

    function getStddev(T_) result(stddev_) ! 設定温度を引数とし、ランダム力計算用の標準偏差を出力する関数
        use parameters
        implicit none
        double precision, intent(in) :: T_
        double precision :: stddev_

        stddev_ = dsqrt(2.000d0 * DAMP * BOLTZ * T_ / dt * 1.0d+33)  ! 無次元 stddev -> 有次元では10^9

    end function getStddev
end module variable

module molecules_struct
    use parameters
    implicit none

    !構造体の定義
    type :: mol_info
        double precision :: pos(3)
        double precision :: vel(3), vtmp(3)
        double precision :: acc(3)
        double precision :: poten, kinet 
    end type mol_info

    type mol_typ
        type(mol_info), allocatable :: mol(:)
    end type mol_typ

    type(mol_typ), dimension(TYPMOL) :: typ  ! 上Pt, 中Ar, 下Pt

end module molecules_struct

! module forPVwin
!     use parameters
!     implicit none
!     integer, parameter :: tlnkoss = nummol(1) + nummol(2) + nummol(3)
!     integer, parameter :: moltype = 1
!     integer, parameter :: ndat = int(stpMax/100)
!     integer, parameter :: ntime0 = 0
!     integer, parameter :: ndt = 1
! end module forPVwin

program main
    use parameters
    use variable
    use molecules_struct
    ! use forPVwin
    implicit none
    double precision :: time
    integer :: i, j

    ! 日付と時刻
    ! character(len=8) :: date
    character(len=100) :: filepath

    filepath = '/home/kawaguchi/' // trim(dir_name)

    ! 配列初期化
    allocate(typ(1)%mol(nummol(1)))
    allocate(typ(2)%mol(nummol(2)))
    allocate(typ(3)%mol(nummol(3)))

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
        open(40,file=trim(filepath) // '/tempe.dat')
        open(41,file=trim(filepath) // '/tempe_PtUp_Layer.dat')
        open(42,file=trim(filepath) // '/tempe_Ar_Layer.dat')
        open(43,file=trim(filepath) // '/tempe_PtLw_Layer.dat')
        
        open(45,file=trim(filepath) // '/tempe_Layer.dat')
    ! 系の周期長さの出力
        open(50,file=trim(filepath) // '/syuuki.dat')
    ! 熱流束のデータ
        open(60,file=trim(filepath) // '/heatflux.dat')
        open(61,file=trim(filepath) // '/pressure.dat')
        
        open(70,file=trim(filepath) // '/force_phantom.dat')
        open(71,file=trim(filepath) // '/force_interface.dat')

    ! 各分子の最終位置データの出力
        open(80,file=trim(filepath) // '/finpos.dat')
    !　分子の色
        ! open(90,file=trim(filepath) // '/mask.dat')

    ! write(15,'(3I7)') moltype, tlnkoss, ndat
    ! do i = 1,ndat
    !     do j = 1, int(nummol(1)/numz(1))
    !         write(90,'(I7)') 15      ! 白色
    !     end do
    !     do j = int(nummol(1)/numz(1)) + 1, nummol(1)
    !         write(90,'(I7)') 14      ! 赤色
    !     end do
    !     do j = 1, nummol(2)
    !         write(90,'(I7)') 7       ! 黄色
    !     end do
    !     do j = 1, int(nummol(3)/numz(3))
    !         write(90,'(I7)') 15      ! 白色
    !     end do
    !     do j = int(nummol(3)/numz(3)) + 1, nummol(3)
    !         write(90,'(I7)') 0       ! 青色
    !     end do
    ! end do
    
    ! condition.datに設計条件を記録
    write(1,*) ''
    write(1, '(A18, F7.4, A3)') 'Scaling Time : ', timeScaling, ' ns '
    write(1, '(A18, F7.4, A3)') 'Relaxation Time : ', timeRelax, ' ns '
    write(1, '(A18, F7.4, A3)') 'Measure Time : ', timeMeasure, ' ns '
    write(1,*) ''
    write(1, '(A8, I5)') 'Ar : ', nummol(2)
    write(1, '(A8, I5, A2)') 'Pt : ', nummol(1), '*2'
    write(1, '(A8, I5)') 'total : ', tlnkoss
    write(1,*) ''
    write(1, '(A23, F4.2)') 'Strength Coefficient : ', angCon
    write(1, '(A23, F5.1)') 'Temperature Ar : ', tempAr
    write(1, '(A23, F5.1)') 'Temperature Pt Upper : ', tempLanPt(1)
    write(1, '(A23, F5.1)') 'Temperature Pt Lower : ', tempLanPt(3)
    write(1,*) ''

    ! ターミナルに表示
    write(6,*) ''
    write(6, '(A18, I8, A4)') 'Scaling Step : ', stpScaling, 'step'
    write(6, '(A18, I8, A4)') 'Relaxation Step : ', stpRelax, 'step'
    write(6, '(A18, I8, A4)') 'Measure Step : ', stpMeasure, 'step'
    write(6,*) ''
    write(6,*) '----------------------- Scaling Step -----------------------'
    write(6, '(I8, A4)') stpScaling, 'step'
    write(6,*) ''
    
    stpNow = 0
    tempLayerAr(:) = 0.0d0
    tempLayerPt(:,:) = 0.0d0

    call seting ! 各分子の初期位置，初期速度などの設定
    call ovito

    do i = 1, stpMax
        stpNow = i
        ! ターミナルに表示
        if(stpNow == stpScaling+1) then
            write(6,*) ''
            write(6,*) '---------------------- Relaxation Step ----------------------'
            write(6,'(I8, A4)') stpRelax, 'step'
            write(6,*) ''
        end if
        if(stpNow == stpScaling + stpRelax+1) then
            write(6,*) ''
            write(6,*) '----------------------- Measure Step ------------------------'
            write(6,'(I8, A4)') stpMeasure, 'step'
            write(6,*) ''

            heatPhantom(:) = 0.0d0 ! パラメータモジュールで初期化するとうまくいかなかったのでここで初期化
            heatInterface(:) = 0.0d0
        end if

        call cpu_time(time)

        ! ステップ数が500の倍数のとき
        if (mod(stpNow,500) == 0) then
            write(6,'(3X, I9, 1X, A3, I10, A1, F15.7, A3)') stpNow, ' / ', stpMax, '(', time, ' s)'
        endif

        ! スケーリング
        if (stpNow <= stpScaling .and. mod(stpNow,100) == 0) then
            call scaling ! 系内の全分子の温度の補正
        endif

        call calc_force_potential ! 各分子に働く力，速度，位置
        if(stpNow > stpScaling + stpRelax) then
            call calc_heatflux
        end if
        call calc_integration

        call bound ! 境界条件の付与

        ! ステップ数が100の倍数+1のとき
        if(mod(stpNow, 100) == 1) then
            ! call record_pos_vel ! 位置，速度を記録
            call record_energy_temp ! エネルギー，温度を記録
            
            if(stpNow <= int((stpScaling + stpRelax/10)/10)) then
                call ovito
            end if
        end if

        if(stpNow > stpScaling + stpRelax) then
            call record_pressure_heatflux ! 熱流束を記録
        end if
    end do

    call calc_force_potential ! 各分子に働く力，速度，位置
    call calc_heatflux
    call calc_integration
    call record_energy_temp
    call record_pressure_heatflux
    call record_final ! 最終状態を記録

    write(6,*) ''
    write(6,*) '----------------------- Completed... ------------------------'
    
    contains

    subroutine seting ! 各分子の初期位置，初期速度などの設定
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: num, i, j, k
        double precision :: coord(3) = 0.0000d0 ! xyz座標
        double precision :: ofst(3)
        double precision :: ran, alpha, beta, cr
        double precision :: v(3)

        do i = 1, COMP
            forCoef(i) = 24.00d0*EPS(i)/SIG(i)  ! 無次元なことに注意 
        end do
            forCoef(4) = angCon*24.00d0*EPS(4)/SIG(4)
        
        syul(:) = syul0(:)
        write(50,*)syul(1)
        write(50,*)syul(2)
        write(50,*)syul(3)
        ! write(15,*)syul0(1), syul0(2), syul0(3)
        ! write(15,*)ntime0, ndt

        do i = 1, TYPMOL
            cutof(:) = syul(:) - CUTOFF*SIG(i)
        end do

        num = 0

        !上段のPt配置
        ofst(1) = STDIST(1)*0.25d0
        ofst(2) = STDIST(1)*0.25d0
        ofst(3) = zsyul0 - STDIST(1)*0.25d0
        do k = 1,numz(1)
            coord(3) = ofst(3) - dble(k-1)*STDIST(1)*0.5d0
            do i = 1,numx(1)
                coord(1) = ofst(1) + dble(i-1)*STDIST(1)*0.5d0
                do j = 1,numy(1)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(1)   !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(1) + STDIST(1)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(1) + STDIST(1)*0.5d0    !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(1)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(1)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        !中段のAr配置
        num = 0
        ofst(1) = syul0(1)*0.50d0 - STDIST(2)*(0.25d0*(numx(2)  -1))
        ofst(2) = syul0(2)*0.50d0 - STDIST(2)*(0.25d0*(numy(2)*2-1))
        ofst(3) = syul0(3)*0.50d0 - STDIST(2)*(0.25d0*(numz(2)  -1))
        do k = 1,numz(2)
            coord(3) = ofst(3) + dble(k-1)*STDIST(2)*0.5d0
            do i = 1,numx(2)
                coord(1) = ofst(1) + dble(i-1)*STDIST(2)*0.5d0
                do j = 1,numy(2)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(2)   !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(2) + STDIST(2)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(2) + STDIST(2)*0.5d0    !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(2)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(2)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        !下段のPt配置
        num = 0
        ofst(1) = STDIST(3)*0.25d0
        ofst(2) = STDIST(3)*0.25d0
        ofst(3) = STDIST(3)*0.25d0
        do k = 1,numz(3)
            coord(3) = ofst(3) + dble(k-1)*STDIST(3)*0.5d0
            do i = 1,numx(3)
                coord(1) = ofst(1) + dble(i-1)*STDIST(3)*0.5d0
                do j = 1,numy(3)
                    if(mod(k,2) == 0) then  !z偶数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(3)   !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(3) + STDIST(3)*0.5d0    !x奇数
                        endif
                    else    !z奇数
                        if(mod(i,2) == 0) then
                            coord(2) = ofst(2) + dble(j-1)*STDIST(3) + STDIST(3)*0.5d0    !x偶数
                        else
                            coord(2) = ofst(2) + dble(j-1)*STDIST(3)   !x奇数
                        endif
                    endif
                    num = num + 1
                    typ(3)%mol(num)%pos(:) = coord(:)
                end do
            end do
        end do

        cr = 1.0d-6
        do j = 1, TYPMOL
            do i = 1, nummol(j)
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

            do i = 1, numx(j)*numy(j)        
                typ(j)%mol(i)%vel(:) = 0.0d0
            end do
        end do
    end subroutine seting

    subroutine scaling ! 系内の全分子の温度の補正
        use parameters
        use variable
        use molecules_struct
        implicit none
        double precision :: temptp, vel2, aimtem, aimnot, baiss
        integer :: i
        integer :: j = 2 ! Arのみ

        temptp = 0.0d0
        do i = 1, nummol(j)
            vel2 = typ(j)%mol(i)%vel(1)**2 + typ(j)%mol(i)%vel(2)**2 + typ(j)%mol(i)%vel(3)**2
            temptp = temptp + vel2
        end do
        temptp = temptp / nummol(j) * 1.0d-16
        aimtem = tempAr
        aimnot = 3.000d0 * BOLTZ * aimtem / MASS(j)
        baiss = dsqrt(aimnot / temptp)

        ! 速度ベクトルのスケーリング
        ! Arのみ
        do i = 1, nummol(j)
            typ(j)%mol(i)%vel(:) = typ(j)%mol(i)%vel(:) * baiss
        end do
    end subroutine scaling

    subroutine ovito
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        double precision :: t
        
        t = dble(stpNow)*dt*1.0d-6 * 100 ! 100ステップおきに計測　単位は [ns]

        write(16, *) tlnkoss
        write(16, '(A8, F4.1, A13, F4.1, A13, F4.1, A39, F6.3)') 'Lattice=', xsyul0, '0.0 0.0 0.0 ', &
                ysyul0, '0.0 0.0 0.0 ', zsyul0, '"Properties="species:S1:pos:R:3" Time=', t
        do j = 1, TYPMOL
            if(j == 2) then
                do i = 1, nummol(j)
                    write(16, '(A2, 3E15.7)') 'Ar', typ(j)%mol(i)%pos(1), typ(j)%mol(i)%pos(2), typ(j)%mol(i)%pos(3)
                end do
            else
                do i = 1, nummol(j)
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
        double precision :: rbook(5000, TYPMOL)
        logical :: judge_same(5000, 5000, TYPMOL), judge_diff(nummol(1), nummol(2), TYPMOL)
        
        double precision :: rnd
        double precision :: div(3), dist
        double precision :: dit2, dit4, dit6, dit8, dit12, dit14
        double precision :: ppp, force, forVec(3)
        double precision :: vel, vmax(5000, TYPMOL), vene(3), sumvene
    
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                typ(j)%mol(i)%acc(:) = 0.0000d0
                typ(j)%mol(i)%poten  = 0.0000d0
                typ(j)%mol(i)%kinet  = 0.0000d0
            end do
        end do
        interForce(:,:,:) = 0.0d0
    
        ! 粒子登録法のために出す速度
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                vel = dsqrt(typ(j)%mol(i)%vel(1)**2 + typ(j)%mol(i)%vel(2)**2 + typ(j)%mol(i)%vel(3)**2)
                if(vmax(i,j) < vel) then
                    vmax(i,j) = vel
                end if
            end do
        end do

        ! 粒子登録法
        ! 最大速度を更新
        if(mod(stpNow,stpUpdate) == 1) then
            do j = 1, TYPMOL
                do i = 1, nummol(j)
                    rbook(i,j) = CUTOFF*SIG(j) + SAFE*stpUpdate*vmax(i,j)*dt
                    vmax(i,j) = 0.0d0
                end do

                do i1 = 1, nummol(j)
                    do i2 = i1+1, nummol(j)
                        div(:) = typ(j)%mol(i1)%pos(:) - typ(j)%mol(i2)%pos(:)
 
                        do k = 1, 2
                            if(div(k) < -syul(k)/2) then
                                div(k) = div(k) + syul(k)
                            else if(div(k) > syul(k)/2) then
                                div(k) = div(k) - syul(k)
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
        ! 同じ分子同士の影響
        do j = 1, TYPMOL
            do i1 = 1, nummol(j)
                do i2 = i1+1, nummol(j)
                    if(.not. judge_same(i1,i2,j)) then
                        cycle
                    else
                        div(:) = typ(j)%mol(i1)%pos(:) - typ(j)%mol(i2)%pos(:)
                        ! カットオフ
                        do k = 1, 3
                            if (div(k) < -cutof(k)) then
                                div(k) = div(k) + syul(k)
                            else if(div(k) > cutof(k)) then
                                div(k) = div(k) - syul(k)
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
                        typ(j)%mol(i1)%poten = typ(j)%mol(i1)%poten + ppp*0.500d0
                        typ(j)%mol(i2)%poten = typ(j)%mol(i2)%poten + ppp*0.500d0
        
                        force  = forCoef(j)*(-2.00d0/dit14+1.0d0/dit8)
                        forVec(:) = -force*div(:)
                        typ(j)%mol(i1)%acc(:) = typ(j)%mol(i1)%acc(:) + forVec(:)/MASS(j)
                        typ(j)%mol(i2)%acc(:) = typ(j)%mol(i2)%acc(:) - forVec(:)/MASS(j)
                    end if
                end do
            end do
        end do
    
        if(mod(stpNow,stpUpdate) == 1) then
            do j = 1, TYPMOL
                if(j == 2) then
                    cycle
                end if
                do i1 = 1, nummol(j) ! Pt
                    do i2 = 1, nummol(2) ! Ar
                        div(:) = typ(j)%mol(i1)%pos(:) - typ(2)%mol(i2)%pos(:)
 
                        do k = 1, 2
                            if(div(k) < -syul(k)/2) then
                                div(k) = div(k) + syul(k)
                            else if(div(k) > syul(k)/2) then
                                div(k) = div(k) - syul(k)
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

        ! 異なる分子同士の影響  ! Ar-Ptの処理    配列は1が上Pt, 2がAr, 3がしたPt
        do j = 1, TYPMOL
            if(j == 2) then
                cycle
            end if
            do i1 = 1, nummol(j)       ! Pt
                do i2 = 1, nummol(2)    ! Ar
                    if(.not. judge_diff(i1,i2,j)) then
                        cycle
                    else
                        div(:) = typ(j)%mol(i1)%pos(:) - typ(2)%mol(i2)%pos(:)
            
                        do k = 1, 3
                            if (div(k) < -cutof(k)) then
                                div(k) = div(k) + syul(k)
                            else if(div(k) > cutof(k)) then
                                div(k) = div(k) - syul(k)
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
                        ppp    = angCon*4.00d0*EPS(4)*(1.0d0/dit12-1.0d0/dit6) ! 異分子間ではangCon(接触角)を忘れずに
                        typ(j)%mol(i1)%poten = typ(j)%mol(i1)%poten + ppp*0.500d0
                        typ(2)%mol(i2)%poten = typ(2)%mol(i2)%poten + ppp*0.500d0
            
                        force  = forCoef(4)*(-2.00d0/dit14+1.0d0/dit8)
                        forVec(:) = -force*div(:)
                        interForce(i1,j,:) = interForce(i1,j,:) - forVec(:) ! 無次元なことに注意　符号が逆な気がする
            
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
    
            do i = numx(j)*numy(j) + 1, 2*numx(j)*numy(j) ! Phantom層のみ
                do k = 1, 3
                    rnd = Random() 
                    ! ランダム力
                    rndForce(i,j,k) = rnd * getStddev(tempLanPt(j)) ! 標準偏差を無次元化したまま扱う
                    ! ダンパー力
                    dmpForce(i,j,k) = -DAMP * typ(j)%mol(i)%vel(k) * 1.0d+14 ! ダンパー力を無次元化したまま扱う
                end do
                ! ランダム力とダンパー力を追加
                typ(j)%mol(i)%acc(:) = typ(j)%mol(i)%acc(:) + (rndForce(i,j,:) + dmpForce(i,j,:)) / MASS(j)*1.0d-3 ! -9+26-20 = -3  (ランダム, ダンパー力の有次元化-9)*(質量の有次元化+26)*(加速度の無次元化-20)
            end do
        end do

        ! 運動エネルギー計算
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                typ(j)%mol(i)%vtmp(:) = typ(j)%mol(i)%vel(:) + typ(j)%mol(i)%acc(:)*0.500d0*dt ! vel(t) = vel(t-dt/2) + acc(t)*dt/2
                sumvene = typ(j)%mol(i)%vtmp(1)**2 + typ(j)%mol(i)%vtmp(2)**2 + typ(j)%mol(i)%vtmp(3)**2
                typ(j)%mol(i)%kinet = 0.500d0*MASS(j)*sumvene
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
    
            do i = numx(j)*numy(j) + 1, 2*numx(j)*numy(j) ! Phantom層  
                heatPhantom(j) = heatPhantom(j) + &
                    &  (rndForce(i,j,1) + dmpForce(i,j,1)) * typ(j)%mol(i)%vtmp(1) + &
                    &  (rndForce(i,j,2) + dmpForce(i,j,2)) * typ(j)%mol(i)%vtmp(2) + &
                    &  (rndForce(i,j,3) + dmpForce(i,j,3)) * typ(j)%mol(i)%vtmp(3) ! 速さの有次元化 10^5
            end do
            
            do i = 1, nummol(j) ! Pt分子全体
                heatInterface(j) = heatInterface(j) + &
                    &  interForce(i,j,1) * typ(j)%mol(i)%vtmp(1) + &
                    &  interForce(i,j,2) * typ(j)%mol(i)%vtmp(2) + &
                    &  interForce(i,j,3) * typ(j)%mol(i)%vtmp(3)
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
                do i = 1, nummol(2)
                    typ(2)%mol(i)%vel(:) = typ(2)%mol(i)%vel(:) + typ(2)%mol(i)%acc(:) * dt   ! vel(t+dt/2) = vel(t-dt/2) + acc(t)*dt
                    typ(2)%mol(i)%pos(:) = typ(2)%mol(i)%pos(:) + typ(2)%mol(i)%vel(:) * dt   ! pos(t+dt)   = pos(t)      + vel(t+dt/2)*dt
                end do
            else
            ! Ptの計算
                ! 固定層
                do i = 1, numx(j)*numy(j)
                    typ(j)%mol(i)%vel(:) = 0.0000d0
                end do
    
                ! その他の層
                do i = numx(j)*numy(j) + 1, int(nummol(j))
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
            do i = 1, nummol(j)
                if(typ(j)%mol(i)%pos(1) < 0.00d0) then
                    typ(j)%mol(i)%pos(1) = typ(j)%mol(i)%pos(1) + syul(1)
                else if(typ(j)%mol(i)%pos(1) > syul(1)) then
                    typ(j)%mol(i)%pos(1) = typ(j)%mol(i)%pos(1) - syul(1)
                endif

                if(typ(j)%mol(i)%pos(2) < 0.00d0) then
                    typ(j)%mol(i)%pos(2) = typ(j)%mol(i)%pos(2) + syul(2)
                else if(typ(j)%mol(i)%pos(2) > syul(2)) then
                    typ(j)%mol(i)%pos(2) = typ(j)%mol(i)%pos(2) - syul(2)
                endif

                ! z方向は周期境界の補正を行わない
                ! if(typ(j)%mol(i)%pos(3) < 0.00d0) then
                !     typ(j)%mol(i)%pos(3) = typ(j)%mol(i)%pos(3) + syul(3)
                ! else if(typ(j)%mol(i)%pos(3) > syul(3)) then
                !     typ(j)%mol(i)%pos(3) = typ(j)%mol(i)%pos(3) - syul(3)
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
    
        do i = 1, nummol(1)
            ! posit_PtUp.dat
            write(10, '(I6, 3E15.7)') (i, typ(1)%mol(i)%pos(k), k = 1,3)
            ! veloc_PtUp.dat
            write(20, '(I6, 3E15.7)') (i, typ(1)%mol(i)%vel(k), k = 1,3)
        end do
        do i = 1, nummol(2)
            ! posit_Ar.dat
            write(11, '(I6, 3E15.7)') (i, typ(2)%mol(i)%pos(k), k = 1,3)
            ! veloc_Ar.dat
            write(21, '(I6, 3E15.7)') (i, typ(2)%mol(i)%vel(k), k = 1,3)
        end do
        do i = 1, nummol(3)
            ! posit_PtLw.dat
            write(12, '(I6, 3E15.7)') (i, typ(3)%mol(i)%pos(k), k = 1,3)
            ! veloc_PtLw.dat
            write(22, '(I6, 3E15.7)') (i, typ(3)%mol(i)%vel(k), k = 1,3)
        end do
    
        do i = numx(1)*numy(1) + 1, 2*numx(1)*numy(1) ! Phantom層
            write(70, '(I6, 6E15.7)') i, rndForce(i,1,1)*1.0d-9, rndForce(i,1,2)*1.0d-9, rndForce(i,1,3)*1.0d-9, &
                                         dmpForce(i,1,1)*1.0d-9, dmpForce(i,1,2)*1.0d-9, dmpForce(i,1,3)*1.0d-9
        end do
    
        do i = 1, nummol(1)!(numz(1)-1)*numx(1)*numy(1) + 1, nummol(1) ! 固液界面層
            write(71, '(I6, 3E15.7)') i, interForce(i,1,1)*1.0d-6, interForce(i,1,2)*1.0d-6, interForce(i,1,3)*1.0d-6 ! 上Pt, Ar, 下Pt
        end do
    
        ! 可視化用
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                ! pos.dat
                write(15, '(3E15.7)') typ(j)%mol(i)%pos(1), typ(j)%mol(i)%pos(2), typ(j)%mol(i)%pos(3)
            end do
        end do
    end subroutine record_pos_vel
    
    subroutine record_energy_temp ! エネルギー，温度を記録
        use parameters
        use variable
        use molecules_struct
        implicit none
        double precision, dimension(TYPMOL) :: totEne, totPot, totKin, temp
        double precision :: allEne, allPot, allKin
        double precision :: kinPtTmp, kinArTmp(numDivAr)
        integer :: cnt(numDivAr)
        double precision :: tempLayerPt_(numz(1),TYPMOL)
        double precision :: tempLayerAr_(numDivAr)
        integer :: i, j, k
        integer :: stp
        stp = stpNow-stpRelax-stpScaling
    
        allEne = 0.0d0
        allPot = 0.0d0
        allKin = 0.0d0
        totEne(:) = 0.0d0
        totPot(:) = 0.0d0
        totKin(:) = 0.0d0
        temp(:) = 0.0d0
        tempLayerPt_(:,:) = 0.0d0
        tempLayerAr_(:) = 0.0d0
    
        do j = 1, TYPMOL
            ! ポテンシャル
            do i = 1, nummol(j)
                totPot(j) = totPot(j) + typ(j)%mol(i)%poten
            end do
            totPot(j) = totPot(j) * 1.0d-16
    
            ! 運動エネルギー
            if (j == 2) then
                ! Ar
                ! do i = 1, nummol(j)
                !     totKin(j) = totKin(j) + typ(j)%mol(i)%kinet
                ! end do
                ! totKin(j) = totKin(j) * 1.0d-16
                ! temp(j) = 2.0d0 * totKin(j) / (3.0d0 * dble(nummol(j)) * BOLTZ)
    
                kinArTmp(:) = 0.0d0
                cnt(:) = 0
                do i = 1, nummol(j)
                    do k = 1, numDivAr
                        if ( (k-1)*zdiv <= (typ(j)%mol(i)%pos(3) - thick(3)) .and. (typ(j)%mol(i)%pos(3) - thick(3)) < k*zdiv) then
                            kinArTmp(k) = kinArTmp(k) + typ(j)%mol(i)%kinet
                            cnt(k) = cnt(k) + 1
                            cycle
                        end if
                    end do
                end do
    
                do k = 1, numDivAr
                    kinArTmp(k) = kinArTmp(k) * 1.0d-16
                    totKin(j) = totKin(j) + kinArtmp(k)
                    tempLayerAr_(k) = 2.0d0 * kinArTmp(k) / (3.0d0 * dble(cnt(k)) * BOLTZ)
                    temp(j) = temp(j) + tempLayerAr_(k)
                end do
    
                temp(j) = temp(j) / dble(numDivAr)
            else
                ! Pt
                do k = 2, numz(j) ! Ptの層の数
                    kinPtTmp = 0.0d0
                    do i = (k-1)*numx(j)*numy(j) + 1, k*numx(j)*numy(j)
                        kinPtTmp = kinPtTmp + typ(j)%mol(i)%kinet
                    end do
    
                    kinPtTmp = kinPtTmp * 1.0d-16
                    totKin(j) = totKin(j) + kinPtTmp
                    tempLayerPt_(k, j) = 2.0d0 * kinPtTmp / (3.0d0 * dble(nummol(j)/numz(j)) * BOLTZ)
                    temp(j) = temp(j) + tempLayerPt_(k, j)
                end do
    
                temp(j) = temp(j) / dble(numz(j)-1)
            end if
    
            totEne(j) = totPot(j) + totKin(j)
            allEne = allEne + totEne(j)
            allPot = allPot + totPot(j)
            allKin = allKin + totKin(j)
        end do
    
        if(stpNow > stpScaling + stpRelax) then
            do j = 1, TYPMOL
                if(j == 2) then
                    do i = 1, numDivAr
                        tempLayerAr(i) = tempLayerAr(i) + tempLayerAr_(i)
                    end do
                else
                    do i = 2, numz(j)
                        tempLayerPt(i,j) = tempLayerPt(i,j) + tempLayerPt_(i,j)
                    end do
                end if
            end do
        end if
        
    
        write(30, '(5E15.7)') stpNow*int(dt)*1.0d-6, totEne(1), totPot(1), totKin(1)  ! energy_PtUp.dat
        write(31, '(5E15.7)') stpNow*int(dt)*1.0d-6, totEne(2), totPot(2), totKin(2)  ! energy_Ar.dat
        write(32, '(5E15.7)') stpNow*int(dt)*1.0d-6, totEne(3), totPot(3), totKin(3)  ! energy_PtLw.dat
        write(35, '(5E15.7)') stpNow*int(dt)*1.0d-6, allEne, allPot, allKin           ! energy_all.dat
        write(40, '(5E15.7)') stpNow*int(dt)*1.0d-6, temp(1), temp(2), temp(3)        ! tempe.dat
    
        ! !!!!!!!!!Pt層を増やすとき必ず変更すること!!!!!!!!!
        ! write(41, '(5E15.7)') stp*int(dt)*1.0d-6, tempLayerPt_(1,1), tempLayerPt_(2,1), tempLayerPt_(3,1), tempLayerPt_(4,1)
        ! write(42, '(16E15.7)') stp*int(dt)*1.0d-6, tempLayerAr_(1), tempLayerAr_(2), tempLayerAr_(3), tempLayerAr_(4)&
        !                                                 , tempLayerAr_(5), tempLayerAr_(6), tempLayerAr_(7), tempLayerAr_(8)&
        !                                                 , tempLayerAr_(9), tempLayerAr_(10), tempLayerAr_(11), tempLayerAr_(12)&
        !                                                 , tempLayerAr_(13), tempLayerAr_(14), tempLayerAr_(15)
        ! write(43, '(5E15.7)') stp*int(dt)*1.0d-6, tempLayerPt_(1,3), tempLayerPt_(2,3), tempLayerPt_(3,3), tempLayerPt_(4,3)
    end subroutine record_energy_temp
    
    subroutine record_pressure_heatflux ! 熱流束を記録
        use parameters
        use variable
        use molecules_struct
        implicit none
        integer :: i, j
        double precision :: heatfluxPhantom(TYPMOL) ! Phantom層からの熱輸送量
        double precision :: heatfluxInterface(TYPMOL) ! 固液界面での熱輸送量
        double precision :: pressure(TYPMOL) ! 圧力
        integer :: stp
        stp = stpNow-stpRelax-stpScaling
    
        heatfluxPhantom(:) = 0.0d0
        heatfluxInterface(:) = 0.0d0
        pressure(:) = 0.0d0
    
        do j = 1, TYPMOL
            if (j == 2) then
                cycle
            end if
    
            heatfluxPhantom(j) = heatPhantom(j) / areaPt * 1.0d+1 ! -9+5+20-15 = +1  (ランダム力-9)*(速度+5)*(面積の逆数+20)*(時間-15)
            heatfluxInterface(j) = heatInterface(j) / areaPt * 1.0d+4 ! -6+5+20-15 = +4  (相互作用力-6)*(速度+5)*(面積の逆数+20)*(時間-15)
            
            do i = 1, nummol(j)
                if(j == 1) then
                    pressure(1) = pressure(1) + interForce(i,1,3) / areaPt * 1.0d+8 ! [MPa]　-6+20-6 = +8  (相互作用力-6)*(面積の逆数+20)*(メガ-6)
                else
                    pressure(3) = pressure(3) - interForce(i,3,3) / areaPt * 1.0d+8
                end if
            end do
        end do
    
        write(60,'(5E15.7)') stp*int(dt)*1.0d-6, heatfluxPhantom(1), heatfluxPhantom(3), heatfluxInterface(1), heatfluxInterface(3)
        write(61,'(4E15.7)') stp*int(dt)*1.0d-6, pressure(1), pressure(2), pressure(3)
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
                do i = 1, numDivAr
                    tempLayerAr(i) = tempLayerAr(i) / dble(stpMeasure/100)
                end do
            else
                do i = 2, numz(j)
                    tempLayerPt(i,j) = tempLayerPt(i,j) / dble(stpMeasure/100)
                end do
            end if
        end do
        
        do j = 1, TYPMOL
            do i = 1, nummol(j)
                ! syuuki.dat
                write(50, '(I6, 6E15.7)') & 
                i, typ(j)%mol(i)%pos(1), typ(j)%mol(i)%pos(2), typ(j)%mol(i)%pos(3), &
                   typ(j)%mol(i)%vel(1), typ(j)%mol(i)%vel(2), typ(j)%mol(i)%vel(3)
            end do
        end do
    
        disz = STDIST(3)*0.25d0
        do i = 2, numz(3)   ! 固定層は除外
            disz = disz + STDIST(3)*0.5d0
            write(45, '(2E15.7)') disz, tempLayerPt(i,3)
        end do
    
        disz = disz - zdiv * 0.5d0
        do i = 1, numDivAr
            disz = disz + zdiv
            write(45, '(2E15.7)') disz, tempLayerAr(i)
        end do
    
        disz = zsyul0 - STDIST(1)*0.25d0
        do i = 2, numz(1)
            disz = disz - STDIST(1)*0.5d0
            write(45, '(2E15.7)') disz, tempLayerPt(i,1)
        end do
    
    end subroutine record_final
end program main