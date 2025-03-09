cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program Local_selfdiffusion_data
! Ar単一バルク系の本計算
! 必要ファイルフォルダ
! 100000
!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 倍精度計算の宣言
      implicit double precision(a-h, o-z)
************************************************************************
! パラメータファイル読み込み
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
! parate.dat  ：パラメータ記入ファイル
! variable.dat：変数記入ファイル
************************************************************************
! 読み込みファイルフォルダ：'in_relax'
      open(1,file='in_restart/finpos_Ar.dat')
      open(2,file='in_restart/finvel_Ar.dat')
! 初期速度のための乱数ファイル読み込み
************************************************************************
************************************************************************
! 出力ファイルフォルダ：'out_relax'
************************************************************************
! 系の周期長さの出力
      open(11,file='out_restart/syuuki.dat')
! 各分子の位置データの出力
      open(12,file='out_restart/posit.dat')
! 各分子の速度データの出力
      open(13,file='out_restart/veloc.dat')
! pvch可視化データの出力
      open(17,file='out_restart/pos.dat')
! 分子の色データの出力
      open(18,file='out_restart/mask.dat')
! 系内のArのエネルギー・温度データの出力
      open(20,file='out_restart/energy_argon.dat')

! msdデータの出力
      open(53,file='out_restart/
     &MeanSquareDisplacement.dat')

! Arの軸方向密度データ
      open(27,file='out_restart/Density_Ar_z.dat')

! 粒子登録法の個数データ記入
      open(36,file='out_restart/book.dat')

! GreenKubo積分データ出力
      open(38,file='out_restart/GreenKubo.dat')

! 台形積分の時間経過の自己拡散係数のデータ出力
      open(39,file='out_restart/SDC_daikei.dat')

! 台形積分用のVACFのデータ出力
      open(41,file='out_restart/VACF.dat')

! 自己拡散係数データ出力
      open(40,file='out_restart/
     &SelfDiffusionCoefficient.dat')

! 並進速度確認データ
      open(60,file='out_restart/veloc_trans.dat')

! 生存確率データ
      open(61,file='out_restart/survive_probability.dat')
! 生存確率を使ったMSDデータ
      open(62,file='out_restart/SP_MSD.dat')
! 生存確率を使ったMSDデータ（確率後処理）
      open(63,file='out_restart/SP_MSD_beta.dat')
! 生存確率を使ったVACFデータ
      open(64,file='out_restart/SP_VACF.dat')
! 生存確率を使ったVACFデータ（確率後処理）
      open(65,file='out_restart/SP_VACF_beta.dat')
! 生存確率を使ったVACFデータ（確率後処理）の台形積分
      open(66,file='out_restart/SP_VACF_beta_daikei.dat')

! 計算終了時の位置・速度データ出力
      open(91,file='out_restart/finpos_Ar.dat')
      open(92,file='out_restart/finvel_Ar.dat')
************************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(6,*)'ReStart Calculation'
      nowstp = 0
      write(36,'(2A10)')'step','Ar-Ar'
************************************************************************
************************************************************************
! 緩和計算終了後の位置・速度データ読み込み
      call read_pos_vel
! 各種パラメータの初期設定
      call seting
      call pvch

!      call record1
!      call record2

! メインルーチン　ここから
      do 100 i_main=1,maxstep
        nowstp = i_main
        if(mod(nowstp,10000) .eq. 0)then
          write(6,*)nowstp
        endif

! 1万ステップ以上で積分前の位置を記録
        if(nowstp .ge. 10000)then
          call positkeep
        endif

! 粒子登録法(40ステップごとに粒子登録)
        if(mod(nowstp,40) .eq. 1)then
          call book
        endif


! 各分子に働く力，速度，位置の分子動力学計算
        call calcu

! 10万ステップ以上で各種拡散係数計算
        if(nowstp .ge. 10000)then
! MSDでの自己拡散係数計算
          call calc_msd
! Green-Kuboでの自己拡散係数計算
! Green-Kuboの速度保存
          call velockeep_GK
! Green-Kuboの内積計算
          call green_kubo
! 生存/非生存の判定
          call survival_hantei
! 生存確率のMSDの計算
          call survival_prob_msd
! 生存確率のVACFの計算
          call survival_prob_vacf
! 生存確率の計算
          call survival_prob

        endif

! 境界条件の付与
        call bound



        if(mod(nowstp, 1000) .eq. 1)then
! データの出力１
!          call record1
! データの出力２
          call record2
!          if(nowstp .le. 100000)then
            call pvch
!          endif
        endif

! 10[ps]ごとに並進速度確認
        if(mod(nowstp, 10000) .eq. 0)then
          call veloc_trans
        endif

! z軸方向密度分布計算
        if(nowstp .ge. 10000)then
          if(mod(nowstp,1000) .eq. 1)then
            call calc_density_dist
          endif
        endif
 100  continue

      call record_msd
      call record_GK
      call record_density
      call record_prob
! 最終位置・速度データ出力
      call record_finish
      write(6,*)'Finish Calculation'
************************************************************************
      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine read_pos_vel
! 粒子位置・速度読み込み
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! Ar分子の緩和計算終了時の位置・速度データ読み込み
      do 101 i1=1, nkoar
        read(1,'(I6,3E24.16)')i,pos1x(i),pos1y(i),pos1z(i)
        read(2,'(I6,3E24.16)')i,vel1x(i),vel1y(i),vel1z(i)
 101  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine seting
! 各種変数・パラメータを初期設定するサブルーチン
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! Ar分子1個の質量
      zmass1 = 1.00D+26*bunsiar/avoga
      sig1   = sigar
      eps1   = epsar
! Ar-Ar
      cforce1 = 24.00D0*eps1/sig1

      xsyul = xsyul0
      ysyul = ysyul0
      zsyul = zsyul0
      write(11,*)xsyul
      write(11,*)ysyul
      write(11,*)zsyul


! Ar-Ar間のカットオフ
      xcutof1 = xsyul - cutoff33*sig1
      ycutof1 = ysyul - cutoff33*sig1
      zcutof1 = zsyul - cutoff33*sig1


      GK_sum1 = 0.000D0
      GK_in_a1bsum=0.000D0
      GK_in_a1tsum=0.000D0
      GK_in_b1sum =0.000D0

! 密度分布計算の初期化
      ndenscount = 0
      do 301 i=1, nmesh
        density_1log(i) = 0.00D0
 301  continue


! msd計算の変位初期化
      do 400 i=1, nshift_m
        do 401 j=1, nkoar
          dsum1x(i,j) = 0.000D0
          dsum1y(i,j) = 0.000D0
          dsum1z(i,j) = 0.000D0
 401     continue
 400  continue

! msd計算の計算ステップカウンタの初期化
      do 402 i=1, nshift_m
        msdcount(i) = 1
        msdcount_st(i) = 1
 402  continue

! msd計算の記録用msd初期化
      do 403 i=1, nmsdtime
        dmsdlog1(i) = 0.000D0
 403  continue

! GreenKubo積分でキープする平均速度の初期化
      do 500 i=1, nGK
        do 501 j=1, nshift_G
          do 502 k1=1, nkoar
            velk_GK1x(i,j,k1) = 0.000D0
            velk_GK1y(i,j,k1) = 0.000D0
            velk_GK1z(i,j,k1) = 0.000D0
 502      continue
 501    continue
 500  continue

! GreenKubo積分での計算ステップカウンタの初期化
      do 503 i=1, nshift_G
        nGKcount(i) = 0
        nGKcount_st(i) = 1
 503  continue

! GreenKubo積分での内積の初期化
      do 504 i=1, nshift_G
        GK_in1(i) = 0.000D0
 504  continue

! GreenKubo積分での記録用変数の初期化
      do 505 i=1, nGK_in
        GK_in1log(i) = 0.000D0
 505  continue

! 生存確率計算用の変数の初期化・設定
! z軸方向のスラブのメッシュ
      z_surv = zsyul/dble(nsurv_mesh_z)
! ステップカウンタ
      do 601 i=1, nshift_s
        nsurvcount(i) = 0
        nsurvcount_st(i) = 1
        do 602 j=1, nkoar
          do 603 k=1, nsurv_mesh_z
            nhantei(i,j,k) = 1
 603      continue
 602    continue
 601  continue
! 個数測定
      do 604 i=1, nshift_s
        do 605 k=1, nsurv_mesh_z
            num_surv(i,k) = 0
            n_slab0(i,k) = 0
 605    continue
 604  continue
! 時間ごとの生存確率保管用配列
      do 606 l=1, nsurv_time
        do 607 k=1, nsurv_mesh_z
            suv_prob(k,l) = 0.0D0
 607    continue
 606  continue

! 生存判定の初期化
      do 608 i=1, nshift_h
        nhancount(i) = 0
        nhancount_st(i) = 1
        do 609 j=1, nkoar
          do 610 k=1, nsurv_mesh_z
            nremain(i,j,k) = 1
 610      continue
 609    continue
 608  continue
      do 611 l=1, nmsdtime
        do 612 m=1, nsurv_mesh_z
          d_svmsdlog1(l,m) = 0.000D0
          d_svvacflog1(l,m) = 0.000D0
 612    continue
 611  continue

      do 613 i=1, nsurv_time
        do 614 k=1, nsurv_mesh_z
            d_svmsdlog2(i,k) = 0.000D0
            d_svvacflog2(i,k) = 0.000D0
 614    continue
 613  continue

! 台形積分用の配列初期化
      do 700 i=1, nGKdai
        GK_dai(i) = 0.00D0
 700  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calcu
! 積分計算のサブルーチン
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! Ar分子に働く力とすべての分子のエネルギーを0に初期化cc
      do 100 i100=1, nkoar
        for1x(i100) = 0.0000D0
        for1y(i100) = 0.0000D0
        for1z(i100) = 0.0000D0
        poten1(i100) = 0.0000D0
        ukine1(i100) = 0.0000D0
 100  continue
************************************************************************

************************************************************************
! Ar-Arの相互作用計算
! Ar分子i1とi2間の各方向の距離を求めるcc
      do 2200 i1= 1, np1
        n1 = np11(i1)
        n2 = np12(i1)
        div1x = pos1x(n1)-pos1x(n2)
        div1y = pos1y(n1)-pos1y(n2)
        div1z = pos1z(n1)-pos1z(n2)
! 周期境界条件の考慮
! もし、分子間の距離がカットオフより大きくても、
! 隣接するセルにある同じ番号の分子との距離を知る必要がある
        if (div1x .lt. -xcutof1) then
          div1x = div1x + xsyul
        elseif(div1x .gt. xcutof1) then
          div1x = div1x - xsyul
        endif

        if (div1y .lt. -ycutof1) then
          div1y = div1y + ysyul
        else if(div1y .gt. ycutof1) then
          div1y = div1y - ysyul
        endif

        if (div1z .lt. -zcutof1) then
          div1z = div1z + zsyul
        elseif(div1z .gt. zcutof1) then
          div1z = div1z - zsyul
        endif

! レナードジョーンズポテンシャルのσで割っておく
! σで割ったカットオフ長さとの比較
! カットオフ長さよりデカければポテンシャルは無視
        div1x = div1x/sig1
        if (div1x .gt. cutoff33) then
          goto 2200
        endif
        if (div1x .lt. -cutoff33) then
          goto 2200
        endif

        div1y = div1y/sig1
        if (div1y .gt. cutoff33) then
          goto 2200
        endif
        if (div1y .lt. -cutoff33) then
          goto 2200
        endif

        div1z = div1z/sig1
        if (div1z .gt. cutoff33) then
          goto 2200
        endif
        if (div1z .lt. -cutoff33) then
          goto 2200
        endif
************************************************************************
! 分子間距離/σの2乗
        dit12 = div1x*div1x + div1y*div1y + div1z*div1z
! 分子間距離
        dist1 = dsqrt(dit12)
! x,y,z方向の距離が全部カットオフより小さくても
! 立体的に距離がカットオフより大きい場合がある
! その時のためのif文
! 先に立体的距離出し手からカットオフするより
! 計算負荷減らせる??
        if(dist1 .gt. cutoff33) then
          goto 2200
        endif
        dit14 = dit12*dit12
        dit16 = dit14*dit12
        dit18 = dit14*dit14
        dit112 = dit16*dit16
        dit114 = dit18*dit16
! pppでレナードジョーンズポテンシャル
        ppp1 = 4.00D0*eps1*(1.00D0/dit112-1.00D0/dit16)
        force1 = cforce1*(-2.00D0/dit114+1.00D0/dit18)
! レナードジョーンズポテンシャルから、各方向の力/質量
! つまり、実質的には加速度の次元で扱う
        force1x = -force1*div1x/zmass1
        force1y = -force1*div1y/zmass1
        force1z = -force1*div1z/zmass1
! ある1対の分子がそれぞれから受ける引力及び斥力
! それらを足し合わせる
        for1x(n1) = for1x(n1) + force1x
        for1x(n2) = for1x(n2) - force1x
        for1y(n1) = for1y(n1) + force1y
        for1y(n2) = for1y(n2) - force1y
        for1z(n1) = for1z(n1) + force1z
        for1z(n2) = for1z(n2) - force1z
! ポテンシャルは2つの分子の相互作用で生じる
! 半分ずつとして簡略化して考える
        poten1(n1) = poten1(n1) + ppp1*0.500D0
        poten1(n2) = poten1(n2) + ppp1*0.500D0
2200  continue
************************************************************************

************************************************************************
! Ar原子の更新
      do 5000 i=1, nkoar
! vneneは各方向の速度線形近似で加速度に時間をかけたものを足す
! 蛙飛び法に従う
        vxene = vel1x(i) + for1x(i)*0.500D0*dt
        vyene = vel1y(i) + for1y(i)*0.500D0*dt
        vzene = vel1z(i) + for1z(i)*0.500D0*dt
! 各方向成分を二乗して足す
        vene = vxene*vxene + vyene*vyene + vzene*vzene
! 運動エネルギーが求められる
        ukine1(i) = 0.500D0*zmass1*vene
5000  continue
! 速度及び位置の積分計算を行う
      do 5100 i=1, nkoar
        vel1x(i) = vel1x(i) + for1x(i)*dt
        vel1y(i) = vel1y(i) + for1y(i)*dt
        vel1z(i) = vel1z(i) + for1z(i)*dt
        pos1x(i) = pos1x(i) + vel1x(i)*dt
        pos1y(i) = pos1y(i) + vel1y(i)*dt
        pos1z(i) = pos1z(i) + vel1z(i)*dt
5100  continue
************************************************************************
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bound
! 周期境界条件付与
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
      do 100 i100=1, nkoar
! x,y軸方向は周期境界条件
        if(pos1x(i100) .lt. 0.00D0)then
          pos1x(i100) = pos1x(i100) + xsyul
        elseif(pos1x(i100) .gt. xsyul)then
          pos1x(i100) = pos1x(i100) - xsyul
        endif
        if(pos1y(i100) .lt. 0.00D0)then
          pos1y(i100) = pos1y(i100) + ysyul
        elseif(pos1y(i100) .gt. ysyul)then
          pos1y(i100) = pos1y(i100) - ysyul
        endif
        if(pos1z(i100) .lt. 0.00D0)then
          pos1z(i100) = pos1z(i100) + zsyul
        elseif(pos1z(i100) .gt. zsyul)then
          pos1z(i100) = pos1z(i100) - zsyul
        endif

 100  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine book
! 粒子登録
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Ar分子の最高速度算出部分
      v1_max = 0.00D0
      do 101 i=1, nkoar
        vel_w=vel1x(i)*vel1x(i)+vel1y(i)*vel1y(i)+vel1z(i)*vel1z(i)
        vel_w=dsqrt(vel_w)
        if(vel_w .gt. v1_max)then
          v1_max=vel_w
        endif
 101  continue

! Ar-Arのbookkeep
      bkr=safe*(cutoff33*sig1+v1_max*dt*40.0D0) !粒子登録する粒子間距離計算
      bkr2=bkr*bkr/sig1/sig1
      bkx = xsyul - bkr
      bky = ysyul - bkr
      bkz = zsyul - bkr
      np1 = 0

! 粒子登録を行うor行わないの判定
      do 201 i=1, nkoar-1
!       周期境界条件を考慮
        do 202 j=i+1, nkoar
          divx = pos1x(i) - pos1x(j)
          if(divx .lt. -xsyul/2.00D0)then
            divx = divx + xsyul
          elseif(divx .gt. xsyul/2.00D0)then
            divx = divx - xsyul
          endif

          divy = pos1y(i) - pos1y(j)
          if(divy .lt. -ysyul/2.00D0)then
            divy = divy + ysyul
          elseif(divy .gt. ysyul/2.00D0)then
            divy = divy - ysyul
          endif

          divz = pos1z(i) - pos1z(j)
          if(divz .lt. -zsyul/2.00D0)then
            divz = divz + zsyul
          elseif(divz .gt. zsyul/2.00D0)then
            divz = divz - zsyul
          endif

!         粒子登録範囲外の粒子をはじく
          divx = divx/sig1
          if(divx .gt. (bkr/sig1))then
            go to 202
          endif
          if(divx .lt. -(bkr/sig1))then
            go to 202
          endif

          divy = divy/sig1
          if(divy .gt. (bkr/sig1))then
            go to 202
          endif
          if(divy .lt. -(bkr/sig1))then
            go to 202
          endif

          divz = divz/sig1
          if(divz .gt. (bkr/sig1))then
            go to 202
          endif
          if(divz .lt. -(bkr/sig1))then
            go to 202
          endif

          dit2 = divx*divx + divy*divy + divz*divz

!           粒子登録を行う粒子を記録
          if(dit2 .lt. bkr2)then
            np1 = np1 + 1
            np11(np1) = i
            np12(np1) = j
          endif
 202    continue
 201  continue


      if(mod(nowstp,1000) .eq. 1)then
        write(36,'(2I10)')nowstp,np1
      endif

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_density_dist
! 軸方向密度分布計算
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! z軸方向の密度のみを測定する
      integer mole1_c(nmesh)
      double precision density_1(nmesh)
************************************************************************
! z軸方向の幅を計算
      z_wall_bot = 0.000D0
      z_wall_top = zsyul
      z_width = (z_wall_top-z_wall_bot)/dble(nmesh)

! 分子の個数を0に初期化
      do 100 i=1, nmesh
        mole1_c(i) = 0
 100  continue

! x,y,z軸方向それぞれについて2種の分子それぞれの個数計算
      do 200 i=1, nmesh
        do 201 i1=1, nkoar
          if(pos1z(i1).ge.((dble(i-1))*z_width+z_wall_bot).and.
     &       pos1z(i1).lt.((dble(i))*z_width+z_wall_bot))then
            mole1_c(i) = mole1_c(i) + 1
          endif
 201    continue
 200  continue

! 個数から密度計算
      do 300 i=1, nmesh
        density_1(i)   = (bunsiar*(dble(mole1_c(i)))/avoga)
     &                   /(xsyul0*ysyul0*z_width*1.0D-30)
        density_1log(i)=density_1log(i) + density_1(i)
 300  continue

! 何回密度分布計算したか
      ndenscount = ndenscount + 1

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine velockeep_GK
! GreenKuboの速度保存
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! GreenKubo積分の計算ステップ
      nGKstp = nowstp - 10000
      nGKstp_st = nGKstp/10000 + 1
      if(nGKstp_st .le. nGK)then
        if(mod(nGKstp,500) .eq. 0)then
          nGKstp_shift = ((nGKstp - (nGKstp_st-1)*10000)/500)+1
          do 101 i1=1, nkoar
            velk_GK1x(nGKstp_st,nGKstp_shift,i1)
     &    = vel1x(i1)
            velk_GK1y(nGKstp_st,nGKstp_shift,i1)
     &    = vel1y(i1)
! z軸方向を除外しておく
!            velk_GK1z(nGKstp_st,nGKstp_shift,i1)
!     &      = vel1z(i1)
 101      continue
        endif
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine green_kubo
! GreenKubo計算
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! 計算ステップ数
      nGKcalcstp = nowstp - 10000

! 内積の計算
      do 200 i=1, nshift_G
        if(nGKcalcstp .ge. (i-1)*500)then
          if(nGKcount(i) .eq. 10000)then
            nGKcount(i) = 0
            nGKcount_st(i) = nGKcount_st(i) + 1
          endif
          if(nGKcount_st(i) .gt. nGK)then
            go to 200
          endif
          GK_in1(i) = 0.000D0
          do 201 i1=1, nkoar
            GK_in1(i) = GK_in1(i)
     &                 +(velk_GK1x(nGKcount_st(i),i,i1)*vel1x(i1)
     &                 + velk_GK1y(nGKcount_st(i),i,i1)*vel1y(i1))
! z軸方向を除外しておく
!     &                   + velk_GK1z(nGKcount_st(i),i,i1)*vel1z(i1))
     &                 /dble(nkoar)
 201      continue
          GK_dai(nGKcount(i)+1) = GK_dai(nGKcount(i)+1)
     &                           +GK_in1(i)/dble(nGK*nshift_G)
          GK_sum1 = GK_sum1 + GK_in1(i)/dble(nGK*nshift_G)
          if(mod(nGKcount(i),100) .eq. 0)then
            nGKcounttime = nGKcount(i)/100 + 1
            GK_in1log(nGKcounttime) = GK_in1log(nGKcounttime)
     &                               +GK_in1(i)/dble(nGK*nshift_G)
          endif
          nGKcount(i) = nGKcount(i) + 1
        endif
 200  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine positkeep
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
      do 600 i1=1, nkoar
        posk1x(i1) = pos1x(i1)
        posk1y(i1) = pos1y(i1)
        posk1z(i1) = pos1z(i1)
 600  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calc_msd
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! 変位計算
      do 101 i1=1, nkoar
        displ1x(i1) = posk1x(i1) - pos1x(i1)
        displ1y(i1) = posk1y(i1) - pos1y(i1)
        displ1z(i1) = posk1z(i1) - pos1z(i1)
 101  continue

! 変位計算ステップ
      ncaldstp = nowstp - 10000

! 変位の和の計算
      do 200 i=1, nshift_m
        if(ncaldstp .ge. (i-1)*500)then
          if(msdcount(i) .eq. 50001)then
            msdcount(i) = 1
            msdcount_st(i) = msdcount_st(i) + 1
            do 301 i1=1, nkoar
              dsum1x(i,i1) = 0.000D0
              dsum1y(i,i1) = 0.000D0
              dsum1z(i,i1) = 0.000D0
 301        continue
          endif
          if(msdcount_st(i) .gt. nmsd)then
            go to 200
          endif
          dmsd1(i) = 0.000D0
          do 201 i1=1, nkoar
            dsum1x(i,i1) = dsum1x(i,i1) + displ1x(i1)
            dsum1y(i,i1) = dsum1y(i,i1) + displ1y(i1)
            dsum1z(i,i1) = dsum1z(i,i1) + displ1z(i1)
            dmsd1(i) = dmsd1(i)
     &                +(dsum1x(i,i1)*dsum1x(i,i1)
     &                + dsum1y(i,i1)*dsum1y(i,i1))/dble(nkoar)
! z軸方向除外
!     &                  + dsum1z(i,i1)*dsum1z(i,i1))/dble(nkoar)
 201      continue

          do 400 j=1, nmsdtime
            if(msdcount(i) .eq. (j*100-50))then
              dmsdlog1(j) = dmsdlog1(j)
     &                     +dmsd1(i)/dble(nmsd*nshift_m)
            endif
 400      continue

          msdcount(i) = msdcount(i) + 1

        endif
 200  continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine veloc_trans
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
      if(nowstp .eq. 10000)then
        write(60,'(4A15)')'time[ns]','vel_x','vel_y','vel_z'
      endif
************************************************************************
      trans_velx = 0.000D0
      trans_vely = 0.000D0
      trans_velz = 0.000D0

      do 101 i1=1, nkoar
        trans_velx = trans_velx + vel1x(i1)
        trans_vely = trans_vely + vel1y(i1)
        trans_velz = trans_velz + vel1z(i1)
 101  continue
! 時間計算[ns]単位
      time_vel_trans = dble(nowstp)/1.00D6

      write(60,'(4D15.7)')time_vel_trans,
     &                    trans_velx,trans_vely,trans_velz

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine survival_prob
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! 計算ステップ数
      nsurvstp = nowstp - 10000
! 残り続けたかの判定
! 0：残り続けなかった(出て行った)
! 1：残り続けた
      do 100 i=1, nshift_s
        if(nsurvstp .ge. (i-1)*500)then
          ! 拡散時間50psの間を計算したとき
          if(nsurvcount(i) .eq. 50000)then
            nsurvcount(i) = 0
            nsurvcount_st(i) = nsurvcount_st(i) + 1
! 一旦全部の粒子が全部の領域に残っていることにする
            do 201 j=1, nkoar
              do 202 k=1, nsurv_mesh_z
                nhantei(i,j,k) = 1
 202          continue
 201        continue
          endif
! 直線的な時間がオーバーしたとき
          if(nsurvcount_st(i) .gt. nsurv)then
            go to 100
          endif
! 残り続けたかの判定エリア
          ! 領域kにある個数を0に初期化
          do 501 k=1, nsurv_mesh_z
            num_surv(i,k) = 0
 501      continue
          ! Ar原子全てについて
          do 301 j=1, nkoar
            ! 各領域において
            do 302 k=1, nsurv_mesh_z
              ! 前のステップでそこにいたら
              if(nhantei(i,j,k) .eq. 1)then
                ! 現在の位置が前のステップでいたところと同じなら
                if((pos1z(j) .ge. dble(k-1)*z_surv) .and.
     &             (pos1z(j) .lt. dble(k)*z_surv))then
                  ! 今その場所にいる個数をプラス1
                  num_surv(i,k) = num_surv(i,k) + 1
                  go to 302
                ! 前のステップでいたところと違うなら
                else
                  ! 判定を0にして出て行ったことにする
                  nhantei(i,j,k) = 0
                endif
              endif
 302        continue
 301      continue
            ! 時刻ゼロでのスラブ内のAr原子数N_slab(0)
          if(nsurvcount(i) .eq. 0)then
              do 401 k=1, nsurv_mesh_z
                n_slab0(i,k) = num_surv(i,k)
 401          continue
          endif

          do 602 k=1, nsurv_mesh_z
            suv_prob(k,(nsurvcount(i)+1)) 
     &    = suv_prob(k,(nsurvcount(i)+1))
     &      + (dble(num_surv(i,k))/dble(n_slab0(i,k)))
     &      / (nshift_s*nsurv) 
 602      continue
 
          nsurvcount(i) = nsurvcount(i) + 1
        endif
        
 100  continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine survival_hantei
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! 計算ステップ数
      nsurvhantei = nowstp - 10000
! 残っているかの判定
! 0：残り続けなかった(出て行った)
! 1：残り続けた
      do 100 i=1, nshift_h
        if(nsurvhantei .ge. (i-1)*500)then
          ! 拡散時間XXpsの間を計算したとき
          if(nhancount(i) .eq. 50000)then
            nhancount(i) = 0
            nhancount_st(i) = nhancount_st(i) + 1
          ! 一旦全粒子が全領域に残っている(存在している)ことにする
            do 201 j=1, nkoar
              do 202 k=1, nsurv_mesh_z
                nremain(i,j,k) = 1
 202          continue
 201        continue
          endif
          ! 直線の時間がオーバーしたとき
          if(nhancount_st(i) .gt. nsurv)then
            go to 100
          endif
          ! 残り続けたか判定
          ! Ar原子すべてについて
          do 301 j=1, nkoar
          ! 各領域において
            do 302 k=1, nsurv_mesh_z
            ! 前のステップでそこにいたら
              if(nremain(i,j,k) .eq. 1)then
              ! 現在の位置が前のステップでいたところと同じなら
                if((pos1z(j) .ge. dble(k-1)*z_surv) .and.
     &             (pos1z(j) .lt. dble(k)*z_surv))then
                  go to 302
                ! 前のステップと違う領域なら
                else
                ! 判定を0にして出て行ったことにする
                  nremain(i,j,k) = 0
                endif
              endif
 302        continue
 301      continue
        endif
      
 100  continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine survival_prob_msd
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
      do 101 i=1, nshift_h      
        if(nsurvhantei .ge. (i-1)*500)then
          if(nhancount(i) .eq. 0)then
           do 201 j=1, nkoar
             do 202 k=1, nsurv_mesh_z
               d_sv_sum1x(i,j,k) = 0.000D0
               d_sv_sum1y(i,j,k) = 0.000D0
               d_sv_sum1z(i,j,k) = 0.000D0
 202         continue
 201       continue
!           do 203 k=1, nsurv_mesh_z
!              n_suv_msd(i,k) = 0
!  203       continue
          endif

          if(nhancount_st(i) .gt. nmsd)then
            go to 101
          endif

          do 301 k=1, nsurv_mesh_z
            n_suv_msd(i,k) = 0
            d_sv_msd(i,k) = 0.000D0
 301      continue

          ! 原子すべてについて
          do 401 j=1, nkoar
            ! 各領域について
            do 402 k=1, nsurv_mesh_z
              if(nremain(i,j,k) .eq. 1)then
                dit_sv1x(i,j,k) = posk1x(j) - pos1x(j)
                dit_sv1y(i,j,k) = posk1y(j) - pos1y(j)
                dit_sv1z(i,j,k) = posk1z(j) - pos1z(j)
                
                d_sv_sum1x(i,j,k) = d_sv_sum1x(i,j,k)
     &                             + dit_sv1x(i,j,k)
                d_sv_sum1y(i,j,k) = d_sv_sum1y(i,j,k)
     &                             + dit_sv1y(i,j,k)
! z軸方向除外
!                d_sv_sum1z(i,j,k) = d_sv_sum1z(i,j,k)
!     &                             + dit_sv1z(i,j,k)
                d_sv_msd(i,k) = d_sv_msd(i,k)
     &                         +d_sv_sum1x(i,j,k)*d_sv_sum1x(i,j,k)
     &                         +d_sv_sum1y(i,j,k)*d_sv_sum1y(i,j,k)
! z軸方向除外
!     &                         +d_sv_sum1z(i,j,k)*d_sv_sum1z(i,j,k)
                n_suv_msd(i,k) = n_suv_msd(i,k) + 1
              endif
 402        continue             
 401      continue
 
          if(nhancount(i) .eq. 0)then
            do 403 k=1, nsurv_mesh_z
                  n_suv0_slab(i,k) = n_suv_msd(i,k)
 403        continue
          endif

          do 501 k=1, nsurv_mesh_z
            ! 残り続けているものの個数で割る
            d_sv_msdA(i,k) = d_sv_msd(i,k)/dble(n_suv_msd(i,k))
            ! N(0)で割る
            d_sv_msdB(i,k) = d_sv_msd(i,k)/dble(n_suv0_slab(i,k))
 501      continue    

          do 601 l=1, nmsdtime
            if(nhancount(i) .eq. (l*100-50))then
              do 602 k=1, nsurv_mesh_z
                d_svmsdlog1(l,k) = d_svmsdlog1(l,k)
     &                            +d_sv_msdA(i,k)/dble(nmsd*nshift_h)
 602          continue
            endif
 601      continue

          do 604 k=1, nsurv_mesh_z
            d_svmsdlog2((nhancount(i)+1),k)
     &    = d_svmsdlog2((nhancount(i)+1),k)
     &     +d_sv_msdB(i,k)/dble(nmsd*nshift_h)
 604      continue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !     nhancount(i) = nhancount(i) + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif
 101  continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine survival_prob_vacf
! 生存確率を使ったVACFの速度保存
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
      do 101 i=1, nshift_h      
        if(nsurvhantei .ge. (i-1)*500)then
          if(nhancount(i) .eq. 0)then
            do 102 j=1, nkoar
              veloc0_1x(i,j) = vel1x(j)
              veloc0_1y(i,j) = vel1y(j)
              veloc0_1z(i,j) = vel1z(j)
 102        continue
          endif

          if(nhancount_st(i) .gt. nmsd)then
            go to 101
          endif

          do 301 k=1, nsurv_mesh_z
            n_suv_vacf(i,k) = 0
            d_sv_vacf(i,k) = 0.000D0
 301      continue
 
        ! 原子すべてについて
          do 401 j=1, nkoar
          ! 各領域について
            do 402 k=1, nsurv_mesh_z
            ! 残り続けているならば
              if(nremain(i,j,k) .eq. 1)then
                d_sv_vacf(i,k) = d_sv_vacf(i,k) 
     &                          +(veloc0_1x(i,j)*vel1x(j)
     &                           +veloc0_1y(i,j)*vel1y(j))
! z軸方向を除外しておく
!     &                           +veloc0_1z(i,j)*vel1z(j))

                n_suv_vacf(i,k) = n_suv_vacf(i,k) + 1
              endif
 402        continue             
 401      continue

          if(nhancount(i) .eq. 0)then
            do 403 k=1, nsurv_mesh_z
              n_suv0_vacf(i,k) = n_suv_vacf(i,k)
 403        continue
          endif

          do 501 k=1, nsurv_mesh_z
            ! 残り続けているものの個数で割る
            d_sv_vacfA(i,k) = d_sv_vacf(i,k)/dble(n_suv_vacf(i,k))
            ! N(0)で割る
            d_sv_vacfB(i,k) = d_sv_vacf(i,k)/dble(n_suv0_vacf(i,k))
 501      continue    

          do 601 l=1, nmsdtime
            if(nhancount(i) .eq. (l*100-50))then
              do 602 k=1, nsurv_mesh_z
                d_svvacflog1(l,k) = d_svvacflog1(l,k)
     &                             +d_sv_vacfA(i,k)/dble(nmsd*nshift_h)
 602          continue
            endif
 601      continue

          do 604 k=1, nsurv_mesh_z
            d_svvacflog2((nhancount(i)+1),k)
     &    = d_svvacflog2((nhancount(i)+1),k)
     &     +d_sv_vacfB(i,k)/dble(nmsd*nshift_h)
 604      continue

          nhancount(i) = nhancount(i) + 1

        endif
 101  continue
      return

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

************************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! 以下、分子運動記述用サブルーチン !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pvch
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
      real    pomx,pomy,pomz
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
      if(nowstp .eq. 0)then
        write(17,'(3I7)') moltype,nmol,ndat1
        write(17,'(3F15.5)') xsyul,ysyul,zsyul
        write(17,'(2I7)') ntime0, ndt
      endif
************************************************************************
      do 500 j=1, nkoar
        pomx = real(pos1x(j))
        pomy = real(pos1y(j))
        pomz = real(pos1z(j))
        write(17,'(3E15.7)') pomx, pomy, pomz
 500  continue

      do 750 i75=1,nkoar
        write(18,'(I7)')nsar
 750  continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine record2
! 粒子のエネルギー・温度記録
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
! 出力ファイルの項目記入
      if(nowstp .eq. 1)then
        write(20,'(5A15)')'time[ps]', 'total_Ar', 'poten_Ar',
     &                    'kine_Ar', 'temp_Ar'
      endif
************************************************************************

! Arのエネルギーの初期化
      toene1  = 0.00D0
      topene1 = 0.00D0
      tokene1 = 0.00D0

*************************************************************************
! Arのエネルギーの合計
      do 101 i1=1, nkoar
        topene1 = topene1 + poten1(i1)
        tokene1 = tokene1 + ukine1(i1)
 101  continue
! Arの総ポテンシャルエネルギー[J]
      topene1 = topene1/1.00D16
! Arの総運動エネルギー[J]
      tokene1 = tokene1/1.00D16
! Arの全エネルギー
      toene1 = topene1 + tokene1
! Arの温度
      temp1 = 2.00D0*tokene1/(3.00D0*dble(nkoar)*boltz)
************************************************************************

! 時間計算[ps]
      time = dble(nowstp)*dt/1000D0
************************************************************************
! Arのエネルギー・温度出力
      write(20,'(5E15.7)')real(time),
     &                    real(toene1), real(topene1),
     &                    real(tokene1), real(temp1)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine record_density
! 密度分布記録
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
      do 101 i=1, nmesh
        density_1log(i)=density_1log(i)/dble(ndenscount)
        zjiku = (dble(i))*z_width - (z_width/2.0D0)
        write(27,'(F15.5,E15.5)')zjiku, real(density_1log(i))
 101  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine record_GK
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
************************************************************************
************************************************************************
      GK_sum_dai = 0.000D0
      do 100 k=1, nsurv_mesh_z
        GK_sum_dai_sp2(k) = 0.000D0
 100  continue
! ファイルに書き込み
      write(38,'(2A15)')'time', 'Ar'
      write(40,'(1A15)')'D_Ar'
! 時間計算
      do 300 i=1, nGK_in
        timeGK = dble(i-1)*0.100D0
        write(38,'(2D15.7)')timeGK,GK_in1log(i)
 300  continue

! 台形積分用のVACF
      do 500 i=1, nGKdai
        time_vacf = dble(i-1)*0.100D0
        write(41,'(2D15.7)')time_vacf,GK_dai(i)
 500  continue

! 台形積分
      do 400 i=1, nGKdai-1
        time_dai = dble(i-1)*0.100D0
        GK_sum_dai = GK_sum_dai + ((GK_dai(i) + GK_dai(i+1))/2.00D0)
     &                           /2.00D0/1.00D5
        write(39,'(2D15.7)')time_dai,GK_sum_dai
 400  continue

! 自己拡散係数のアンサンブル平均出力
! 単位は[m^2/s]
! xy平面のときは2で割る
      sdc1 = (GK_sum1/2.00D0)/1.00D5
      write(40,'(1E15.7)')real(sdc1)

! 生存確率を使ったVACFのデータ出力
      do 201 i=1, nmsdtime
! vacfの時間計算（ピコ秒単位）
        time_suvvacf = dble(i-1)/10.0D0
      write(64,'(9E15.7)')real(time_suvvacf), 
     &                      real(d_svvacflog1(i,1)),
     &                      real(d_svvacflog1(i,2)),
     &                      real(d_svvacflog1(i,3)),
     &                      real(d_svvacflog1(i,4)),
     &                      real(d_svvacflog1(i,5)),
     &                      real(d_svvacflog1(i,6)),
     &                      real(d_svvacflog1(i,7)),
     &                      real(d_svvacflog1(i,8))
 201  continue
! 生存確率を使ったVACFのデータ出力（確率で後処理型）
      do 202 i=1, nsurv_time
!         VACFの時間計算(ピコ秒単位)
        time_suvvacf = dble(i-1)/1000.0D0
        do 203 k=1, nsurv_mesh_z
          dvacf_pk(i,k) = d_svvacflog2(i,k) 
     &                  / suv_prob(k, i)
 203    continue 
        write(65,'(9E15.7)')real(time_suvmsd), 
     &                      real(dvacf_pk(i,1)),real(dvacf_pk(i,2)),
     &                      real(dvacf_pk(i,3)),real(dvacf_pk(i,4)),
     &                      real(dvacf_pk(i,5)),real(dvacf_pk(i,6)),
     &                      real(dvacf_pk(i,7)),real(dvacf_pk(i,8))
 202  continue

! 台形積分
      do 401 i=1, nsurv_time-1
        time_dai = dble(i-1)*0.00100D0
        do 402 k=1, nsurv_mesh_z
          GK_sum_dai_sp2(k) = GK_sum_dai_sp2(k) 
     &   +((d_svvacflog2(i,k) + d_svvacflog2(i+1,k))/2.00D0)
     &                           /2.00D0/1.00D5
 402    continue
        write(66,'(9E15.7)')real(time_dai),
     &                      real(GK_sum_dai_sp2(1)),
     &                      real(GK_sum_dai_sp2(2)),
     &                      real(GK_sum_dai_sp2(3)),
     &                      real(GK_sum_dai_sp2(4)),
     &                      real(GK_sum_dai_sp2(5)),
     &                      real(GK_sum_dai_sp2(6)),
     &                      real(GK_sum_dai_sp2(7)),
     &                      real(GK_sum_dai_sp2(8))
 401  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine record_msd
! MSDの記録
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
************************************************************************
      double precision time_msd(nmsdtime)
************************************************************************
! ファイルに書き込み
      write(53,'(2A15)')'time[ps]', 'Ar_MSD'


! msdのデータ出力
      do 200 i=1, nmsdtime
!         msdの時間計算(ピコ秒単位)
        time_msd(i) = dble(i)/10.0D0 - 0.0500D0
        write(53,'(2D15.7)')time_msd(i), dmsdlog1(i)
 200  continue

! 生存確率を使ったMSDのデータ出力
      do 201 i=1, nmsdtime
!         msdの時間計算(ピコ秒単位)
        time_suvmsd = dble(i)/10.0D0 - 0.0500D0
        write(62,'(9E15.7)')real(time_suvmsd), 
     &                      real(d_svmsdlog1(i,1)),
     &                      real(d_svmsdlog1(i,2)),
     &                      real(d_svmsdlog1(i,3)),
     &                      real(d_svmsdlog1(i,4)),
     &                      real(d_svmsdlog1(i,5)),
     &                      real(d_svmsdlog1(i,6)),
     &                      real(d_svmsdlog1(i,7)),
     &                      real(d_svmsdlog1(i,8))
 201  continue

! 生存確率を使ったMSDのデータ出力（確率で後処理型）
      do 202 i=1, nsurv_time
!         msdの時間計算(ピコ秒単位)
        time_suvmsd = dble(i)/1000.0D0 - 0.000500D0
        do 203 k=1, nsurv_mesh_z
          dmsd_pk(i,k) = d_svmsdlog2(i,k) 
     &                  / suv_prob(k, i)
 203    continue 
        write(63,'(9E15.7)')real(time_suvmsd), 
     &                      real(dmsd_pk(i,1)),real(dmsd_pk(i,2)),
     &                      real(dmsd_pk(i,3)),real(dmsd_pk(i,4)),
     &                      real(dmsd_pk(i,5)),real(dmsd_pk(i,6)),
     &                      real(dmsd_pk(i,7)),real(dmsd_pk(i,8))
 202  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine record_prob
! 生存確率の記録
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
************************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(61,'(9A15)')'time[ps]',
     &             'region1','region2','region3','region4',
     &             'region5','region6','region7','region8'
      do 101 l=1, nsurv_time
        survtime = dble(l-1)*1.0D-3
        write(61,'(9E15.7)')REAL(survtime),
     &         real(suv_prob(1,l)), real(suv_prob(2,l)),
     &         real(suv_prob(3,l)), real(suv_prob(4,l)),
     &         real(suv_prob(5,l)), real(suv_prob(6,l)),
     &         real(suv_prob(7,l)), real(suv_prob(8,l))
 101  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine record_finish
! 最終位置・速度の記録
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      include 'in_restart/parate.dat'
      include 'in_restart/variable.dat'
************************************************************************
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 101 i=1, nkoar
        write(91,'(I6,3E24.16)')i,pos1x(i),pos1y(i),pos1z(i)
        write(92,'(I6,3E24.16)')i,vel1x(i),vel1y(i),vel1z(i)
 101  continue

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc