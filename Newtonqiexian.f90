program newton
	real*8 a,b,t,dt
	a=0.0;b=3.0
	t=1.0;e=0.01
	dt=1.0
	
	open(13,file="Newton")
	write(13,"(A12,A12)")"t-t0","t"
	do while(dt>=e)
		dt=(((3*t**2-2)/(6*t))**2)**0.5
		t=t-(3*t**2-2)/(6*t)
		write(13,"(f12.6,f12.6)")dt,t
	end do
	write(13,*)t,t**3-2*t+1,e
end program newton
