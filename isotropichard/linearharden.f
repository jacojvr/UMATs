      subroutine fisotropic(sy,dsy,deqpl,dtime,temp,
     & statev,tempstatev,nstatv,props,nprops)
      implicit real*8(a-h,o-z)
      dimension statev(nstatv),tempstatev(nstatv),props(nprops)
      eqpl = statev(1)+deqpl
      sy = props(3)+props(4)*eqpl
      dsy = props(4)
      tempstatev(1) = eqpl
      return
      end