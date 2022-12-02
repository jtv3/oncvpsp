subroutine writeorb( prefix, zz, nn, ll, rr, uu, irc )
 implicit none
 integer, parameter :: dp=kind(1.0d0)

  character(len=4), intent( in )  :: prefix
  real(DP), intent( in ) :: zz
  integer, intent( in ) :: nn, ll, irc
  real(DP) :: rr(1:irc), uu(1:irc)

  integer :: ii
  character( len=20) :: fnam

   write( fnam,  '(2a4,1i3.3,2(1a1,1i2.2))' ) prefix, 'orbz', nint( zz ), 'n', nn, 'l', ll
   write( 6, * ) fnam, irc, sizeof( rr ), sizeof(uu)

   open( unit=99, file=fnam, form='formatted', status='unknown' )
   rewind 99
   write ( 99, '(1a1,1i8)' ) '#', irc
   do ii = 1, irc
      write ( 99, '(2(1x,1e22.15))' ) rr( ii ), uu( ii )
   end do
   close( unit=99 )

end subroutine writeorb
