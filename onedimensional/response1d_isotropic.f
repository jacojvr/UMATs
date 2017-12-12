      subroutine risotropic(stress, dstrain, dtime, temp,
     & statev, nstatv, props, nprops)
!f2py real*8 intent(in,out) :: stress,backstress
!f2py real*8 intent(in) :: dstrain,dtime,temp
!f2py integer intent(hide),depend(statev) :: nstatv = len(statev)
!f2py real*8 intent(in,out),dimension(nstatv) :: statev
!f2py integer intent(hide),depend(props) :: nprops = len(props)
!f2py real*8 intent(in),dimension(nprops) :: props
!
      implicit real*8(a-h,o-z)
      dimension statev(nstatv),props(nprops),tempstatev(nstatv)
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, six=6.d0,
     & newtonmax=10, toler=1.0d-8)
!
! elastic properties
      call shearmod(eg,enu,temp,props,nprops)
      eg2=two*eg
      emod = eg2*(one+enu)
!
! initialize temporary statevariable array
      do i =1,nstatv
        tempstatev(i) = statev(i)
      end do
!
! calculate predictor stress
      stress=stress+emod*dstrain
      strial=dabs(stress)
!
! check yield condition
      call fisotropic(sy,dsy,zero,dtime,temp,
     & statev,tempstatev,nstatv,props,nprops)
!
! determine if actively yielding
      residual = strial-sy
      if(dabs(residual).gt.toler)then
!
! direction of the plastic flow
      flow=stress/strial
!
! newton raphson return mapping
      cop=toler
      kount = 1
      do while((dabs(residual).gt.toler).and.(kount.lt.newtonmax))
        cop=dabs(cop+residual/(emod+dsy))
        deqpl = cop
        call fisotropic(sy,dsy,deqpl,dtime,temp,
     & statev,tempstatev,nstatv,props,nprops)
        residual=strial-emod*cop-sy
      end do
!
! update stress solution
      stress=flow*sy
      endif
!
! update state variables
      do k1=1,nstatv
        statev(k1) = tempstatev(k1)
      end do
      return
      end
!