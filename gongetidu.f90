program gonge
	use gongesub1
	real*8 x(2),A(2,2),b(2)
	integer n
	write(*,"(6A13)")"f","x1","x2","g1","g2","i"
	data a /4.0d0,-1.0d0,-1.0d0,2.0d0/
	data b/0.0d0,0.0d0/
	data x/1.0d0,1.0d0/
	n=2
	call cg(A,b,x,n)
end program gonge
