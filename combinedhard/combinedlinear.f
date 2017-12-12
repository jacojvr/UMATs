      subroutine fcombined(sy,dsy,sb,dsb,deqpl,dtime,temp,
     & statev,tempstatev,nstatv,props,nprops)
      implicit real*8(a-h,o-z)
      dimension statev(nstatv),tempstatev(nstatv),props(nprops)
      eqpl = statev(7)+deqpl
      sy = props(3)+props(4)*eqpl
      dsy = props(4)
      sb = props(5)*eqpl
      dsb = props(5)
      tempstatev(7) = eqpl
      return
      end