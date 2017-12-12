      subroutine cyclreversal(statev,nstatv)
      implicit real*8(a-h,o-z)
      dimension statev(nstatv)
      ! density - recoverable
      statev(2) = statev(2)-statev(3)
      ! reset recoverable
      statev(3) = 0.d0
      ! swap backstress
      statev(4) = -statev(4)
      return
      end