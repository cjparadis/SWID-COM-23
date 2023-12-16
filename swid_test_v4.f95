program swid_com_23

	implicit none

	! Stress period 1: forced-gradient radial transport
	! Governing differential equation: R*dC/dt=Dr*d^2C/dr^2-vr*dC/dr+kC
	! where: Dr=alphar*vr
	! Disc. adv. term: vr*dC/dr=(v[j]+v[j-1])/2*(C[j,n]-C[j-1,n])/delr
	! where: v[j]=Q/(2*pi*r[j]*b*n)
	! Disc. disp. term: alphar*(((v[j+1]+v[j])/2*(C[j+1,n]-C[j,n])/delr)-((v[j]+v[j-1])/2*(C[j,n]-C[j-1,n])/delr))/delr
	! Disc. temp. term: R*dC/dt=R*(C[j,n+1]-C[j,n])/delt
	! Disc. Rxn. term: k*C[j,n]
	
	! Stress period 2: natural-gradient horizontal transport
	! Governing differential equation: RdC/dt=Dx*d^2C/dx^2-vx*dC/dx+kC 
	! where: Dx=alphax*vx
	! Disc. adv. term: vx*dC/dx=vx*(C[j,n]-C[j-1,n])/delx
	! where: vx=-Ki/n	
	! Disc. disp. term: alphax*vx(C[j+1,n]-2*C[j,n]+C[j-1,n])/delx^2
	! Disc. temp. term: R*dC/dt=R*(C[j,n+1]-C[j,n])/delt
	! Disc. Rxn. term: k*C[j,n]
	
	! Let alphar=alphax=alpha & delr=delx
	! Units: feet, days
	! Model specs: 
		! 1) One dimensional in radius (stress period 1) and horizontal (stress period 2)
		! 2) 50 feet in length
		! 3) 1-inch grid space or 1/12-foot grid space
		! 4) 10^-5 days time space (stress period 1), 10^-3 days time space (stress period 2)
		! 5) Specificed concentration boundary conditions at inlet and outlet
		! 6) Stress period 1: inlet conc. = injection conc., outlet conc. = aquifer conc., initial conc. = uniform aquifer conc.
		! 7) Stress period 2: inlet and outlet concs. = aquifer conc., initial conc. = varibale aquifer conc., i.e., Gaussian distribution
		! 8) Transport equations solved by explicit finite difference, similar to MT3DMS Group C solutions in Table 1 on pg. 36 (1999)
	
	! Declare, Define, Compute, Output

	! Declare
	integer :: i, j, n, nend, kount1, kount2, kount3, kount4, eye, aye
	integer,parameter :: col=301, injw=301, sim=3
	real,parameter :: pi=3.14159, e=2.178
	
	! Stress period 1 (S1)
	real,dimension(sim,col) :: ColdS1, CnewS1 ! 3 rows or sims & 301 columns or 25 feet (1/12ft*301-1=25ft)
	real,dimension(col) :: r, vr ! 301 columns or 25 feet
	real,dimension(col,sim+1) :: Cr ! 301 rows or 25 feet & 4 columns of distance plus 3 sims
	real,parameter :: a1=0.278393, a2=0.230389, a3=0.000972, a4=0.078108 ! Constants for erfc for analytical solution
	character(len=20),dimension(sim+1+sim) :: HeaderS1 ! 7 columns, time, 3 numerical sims, and 3 analytical sims
	character(len=3) :: AM ! Analytical model, let user decide on this...for now 
	
	! Stress period 2 (S2)
	real,dimension(sim,col+col-1) :: ColdS2, CnewS2 ! 3 rows or sims & 601 columns or 50 feet (1/12ft*601-1=50ft)
	real,dimension(col+col-1) :: x ! 601 columns or 50 feet
	real,dimension(col+col-1,sim+1+sim) :: Cx ! 601 rows or 50 feet with 7 columns, 1 distance, 3 numerical sims,and 3 analytical sims
	character(len=20),dimension(sim+1) :: HeaderS2 ! 4 columns of time plus 3 sims
	
	! Breakthrough curve
	real,dimension(col-1,sim+1) :: Ct ! 300 rows of time with 4 columns of time plus 3 sims
	real :: Qi, b, theta, d, ti, Ci, Ca ! User inputs/knowns
	real,dimension(sim) :: vx, alpha, Ret, k ! 3 velocity, 3 alpha, 3 retention, & 3 decay simulations (sims)
	real :: delr, delx, deltS1, deltS2, time, vol
	real :: Pe1, Pe2, Pe3, rmin ! Stability parameters 
	real :: AdvStab1, AdvStab2, AdvStab3, AdvStab4 ! Stability parameters
	real :: DispStab1, DispStab2, DispStab3, DispStab4, DispStab5, DispStab6 ! Stability parameters
	real :: DecayStab1, DecayStab2, DecayStab3 ! Stability parameters
	
	! Mass budget
	real,dimension(3) :: Mo, Mt, So, Si, delM, Disc ! Mass old, new, source, sink, change, and discrepancy
	real,dimension(100,sim+1+sim) :: DiscS1, DiscS2 ! Detailed mass discrepancy for stress periods 1 and 2
	character(len=20),dimension(sim+1+sim) :: HeaderDiscS1, HeaderDiscS2 ! 7 columns of time plus disc and delM sims
	character(len=3) :: MB ! Let user decide on this...for now

	! Other stuff
	real :: A, AA, AAA, sumA, sumAA ! Arithmetic placeholders

	! Default input parameters from Hoss (2022) thesis (https://dc.uwm.edu/etd/2900/)
	Qi=100./(7.5*60.) !gpm
	ti=7.5 !hours
	b=5. !ft
	theta=0.38 !porosity
	Ci=500. !mg/L
	Ca=0.5 !mg/L
	
	! Three simulations
	vx(1)=0.5
	vx(2)=1.
	vx(3)=2.
	
	do i=1,sim
		alpha(i)=0.1
		Ret(i)=1.
		k(i)=0.
	end do
	
	! Analytical model and mass budget
	AM="Yes"
	MB="Yes"
	
	! Header for program
	write(*,*) ""
	write(*,*) "============SWID-COM-23============"
	write(*,*) "------------NSF Uranium Grant------"
	write(*,*) "************Award #2229869*********"
	write(*,*) "............UW-Milwaukee..........."
	write(*,*) ""
	
	! Get basic input parameters (Stress Period 1)
	write(*,*) "ENTER INPUT PARAMETERS FOR INJECTION PHASE:"
	write(*,*) ""
	write(*,*) "Pumping rate (gpm)="
	read(*,*) Qi
	write(*,*) "Pumping time (hours)="
	read(*,*) ti
	write(*,*) "Aquifer thickness (ft)="
	read(*,*) b
	write(*,*) "Aquifer porosity (-)="
	read(*,*) theta
	write(*,*) "Solute concentration injection fluid (mg/L)="
	read(*,*) Ci
	write(*,*) "Solute concentration aquifer fluid (mg/L)="
	read(*,*) Ca

	! Get basic input parameters (Stress Period 2)
	write(*,*) ""
	write(*,*) "ENTER INPUT PARAMETERS FOR DRIFT PHASE:"
	write(*,*) ""
	write(*,*) "Groundwater velocities for simulations 1, 2, & 3 (ft/day)="
	read(*,*) vx(1),vx(2),vx(3)
	write(*,*) "Aquifer dispersivities for simulations 1, 2, & 3 (ft)="
	read(*,*) alpha(1), alpha(2), alpha(3)
	write(*,*) "Solute retentions for simulations 1, 2, & 3 (-)="
	read(*,*) Ret(1), Ret(2), Ret(3)
	write(*,*) "Solute decay rates for simulations 1, 2, & 3 (1/day)="
	read(*,*) k(1), k(2), k(3)
	write(*,*) ""
	
	! Compare injection phase to analytical
	write(*,*) ""
	write(*,*) "ANALYTICAL MODEL:"
	write(*,*) ""
	write(*,*) "Enter Yes or No for analytical model"
	read(*,*) AM
	
	! Mass budget yes or no
	write(*,*) ""
	write(*,*) "MASS BUDGET:"
	write(*,*) ""
	write(*,*) "Enter Yes or No for mass budget"
	read(*,*) MB
	
	! k is a decay rate, must ALWAYS be negative
	do i=1,sim
		k(i)=-abs(k(i))
	end do
	
	! Define
	Qi=Qi*192.5 ! gpm to ft^3/day
	ti=ti*0.04166667 ! hours to days
	
	d=1./12. ! Grid spacing of one inch
	
	delr=d ! Grid space
	delx=d ! Grid space
	
	deltS1=0.00001 ! Time step stress period 1 (injection)
	deltS2=0.001 ! Times step stress period 2 (drift)
	
	! Conditionally stable parameters
	Pe1=((Qi/(2.*pi*delr*b*theta))*delr)/((Qi/(2.*pi*delr*b*theta))*alpha(1)) ! Peclet needs to be less than 4 for stability
	Pe2=((Qi/(2.*pi*delr*b*theta))*delr)/((Qi/(2.*pi*delr*b*theta))*alpha(2)) ! Peclet needs to be less than 4 for stability
	Pe3=((Qi/(2.*pi*delr*b*theta))*delr)/((Qi/(2.*pi*delr*b*theta))*alpha(3)) ! Peclet needs to be less than 4 for stability
	
	rmin=((Qi*ti)/(pi*b*theta))**(1./2.) ! Radius of influence needs to be less than ~10 feet for no boundary effects
	
	AdvStab1=Ret(1)/((Qi/(2.*pi*delr*b*theta))/delr) ! See MT3DMS Guide Pg. 54
	
	AdvStab2=Ret(1)/(vx(1)/delx) ! See MT3DMS Guide Pg. 54
	AdvStab3=Ret(2)/(vx(2)/delx) ! See MT3DMS Guide Pg. 54
	AdvStab4=Ret(3)/(vx(3)/delx) ! See MT3DMS Guide Pg. 54
	
	DispStab1=(0.5*Ret(1))/((alpha(1)*(Qi/(2.*pi*delr*b*theta)))/(delr**2.)) ! See MT3DMS Guide Pg. 54
	DispStab2=(0.5*Ret(2))/((alpha(2)*(Qi/(2.*pi*delr*b*theta)))/(delr**2.)) ! See MT3DMS Guide Pg. 54
	DispStab3=(0.5*Ret(3))/((alpha(3)*(Qi/(2.*pi*delr*b*theta)))/(delr**2.)) ! See MT3DMS Guide Pg. 54
	
	DispStab4=(0.5*Ret(1))/((alpha(1)*vx(1))/(delx**2.)) ! See MT3DMS Guide Pg. 54
	DispStab5=(0.5*Ret(2))/((alpha(2)*vx(2))/(delx**2.)) ! See MT3DMS Guide Pg. 54
	DispStab6=(0.5*Ret(3))/((alpha(3)*vx(3))/(delx**2.)) ! See MT3DMS Guide Pg. 54
	
	DecayStab1=1./abs(k(1)) ! See MT3DMS Guide Pg. 54
	DecayStab2=1./abs(k(2)) ! See MT3DMS Guide Pg. 54
	DecayStab3=1./abs(k(3)) ! See MT3DMS Guide Pg. 54
	
	!Decrease time steps during Injection Phase or Drift Phase if not stable (turned off for now...)
	
	!if(deltS1>AdvStab1) then
	!		deltS1=deltS1/(10.*deltS1/AdvStab1)
	!	else if(deltS1>DispStab2) then
	!		deltS1=deltS1/(10.*deltS1/DispStab1)
	!	else if(deltS1>DispStab3) then
	!		deltS1=deltS1/(10.*deltS1/DispStab2)
	!	else if(deltS1>DispStab4) then	
	!		deltS1=deltS1/(10.*deltS1/DispStab3)
	!end if
	
	! WARNINGS START :( ======================================================================== :(
	if(Pe1>=4) then
		write(*,*) ""
		write(*,*) " WARNING #1: Peclet # >= 4 in simulation 1, numerical dispersion/artificial oscillation"
		write(*,*) ""
	else if(Pe2>=4) then
		write(*,*) ""
		write(*,*) " WARNING #1: Peclet # >= 4 in simulation 2, numerical dispersion/artificial oscillation"
		write(*,*) ""
	else if(Pe3>=4) then
		write(*,*) ""
		write(*,*) " WARNING #1: Peclet # >= 4 in simulation 3, numerical dispersion/artificial oscillation"
		write(*,*) ""
	end if
	
	if(rmin>=10.) then
		write(*,*) ""	
		write(*,*) " WARNING #2: Minimum radius of influence ~>= 10 feet, boundary condition effect"
		write(*,*) ""
	end if

	if(rmin<0.5) then
		write(*,*) ""
		write(*,*) " WARNING #3: Minimum radius of influence ~< 0.5 feet, boundary condition effect"
		write(*,*) ""
	end if

	if(deltS1>AdvStab1) then
		write(*,*) ""
		write(*,*) " WARNING #4: Advection term unstable for injection phase"
		write(*,*) ""
	end if

	if(deltS2>AdvStab2)	then
		write(*,*) ""
		write(*,*) " WARNING #5: Advection term unstable for simulation 1 drift phase"
		write(*,*) ""
	else if(deltS2>AdvStab3) then
		write(*,*) ""
		write(*,*) " WARNING #5: Advection term unstable for simulation 2 drift phase"
		write(*,*) ""	
	else if(deltS2>AdvStab4) then
		write(*,*) ""
		write(*,*) " WARNING #5: Advection term unstable for simulation 3 drift phase"
		write(*,*) ""
	end if
	
	if(deltS1>DispStab1) then
		write(*,*) ""	
		write(*,*) " WARNING #6: Dispersion term unstable for simulation 1 injection phase"
		write(*,*) ""
	else if(deltS1>DispStab2) then	
		write(*,*) ""	
		write(*,*) " WARNING #6: Dispersion term unstable for simulation 2 injection phase"
		write(*,*) ""
	else if(deltS1>DispStab3) then	
		write(*,*) ""	
		write(*,*) " WARNING #6: Dispersion term unstable for simulation 3 injection phase"
		write(*,*) ""
	end if
	
	if(deltS2>DispStab4) then
		write(*,*) ""	
		write(*,*) " WARNING #7: Dispersion term unstable for simulation 1 drift phase"
		write(*,*) ""
	else if (deltS2>DispStab5) then	
		write(*,*) ""	
		write(*,*) " WARNING #7: Dispersion term unstable for simulation 2 drift phase"
		write(*,*) ""
	else if (deltS2>DispStab6) then
		write(*,*) ""	
		write(*,*) " WARNING #7: Dispersion term unstable for simulation 3 drift phase"
		write(*,*) ""
	end if
	
	if(deltS1>DecayStab1) then	
		write(*,*) ""	
		write(*,*) " WARNING #8: Decay term unstable for simulation 1 injection phase"
		write(*,*) ""
	else if(deltS1>DecayStab2) then	
		write(*,*) ""	
		write(*,*) " WARNING #8: Decay term unstable for simulation 2 injection phase"
		write(*,*) ""
	else if(deltS1>DecayStab3) then	
		write(*,*) ""	
		write(*,*) " WARNING #8: Decay term unstable for simulation 3 injection phase"
		write(*,*) ""
	end if
	
	if(deltS2>DecayStab1) then	
		write(*,*) ""	
		write(*,*) " WARNING #9: Decay term unstable for simulation 1 drift phase"
		write(*,*) ""
	else if(deltS2>DecayStab2) then	
		write(*,*) ""	
		write(*,*) " WARNING #9: Decay term unstable for simulation 2 drift phase"
		write(*,*) ""
	else if(deltS2>DecayStab3) then	
		write(*,*) ""	
		write(*,*) " WARNING #9: Decay term unstable for simulation 3 drift phase"
		write(*,*) ""
	end if
	! WARNINGS END :) ========================================================================== :)
	
	! Fill arrays with zero place holders
	do i=1,sim
		do j=1,col
			ColdS1(i,j)=0.
			CnewS1(i,j)=0.
			r(j)=0.
		end do
	end do
	
	do j=1,sim+1
		do i=1,col
			Cr(i,j)=0.
		end do
	end do
	
	do i=1,sim
		do j=1,col+col-1
			ColdS2(i,j)=0.
			CnewS2(i,j)=0.
			x(j)=0.
		end do
	end do
	
	do j=1,sim+1+sim
		do i=1,col+col-1
			Cx(i,j)=0.
		end do
	end do
	
	do j=1,sim+1
		do i=1,col-1
			Ct(i,j)=0.
		end do
	end do
	
	do j=1,sim+1
		do i=1,100
			DiscS1(i,j)=0.
		end do
	end do		
	
	! Set boundary conditions (stress period 1) at inlet (left side) and outlet (right side)
	do i=1,sim
		ColdS1(i,1)=Ci ! Inlet
		CnewS1(i,1)=Ci ! Inlet
	end do	
	
	do i=1,sim
		ColdS1(i,col)=Ca ! Outlet
		CnewS1(i,col)=Ca ! Outlet
	end do

	! Set initial condition (stres period 1) within boundaries
	do i=1,sim
		do j=2,col-1
			ColdS1(i,j)=Ca
		end do
	end do
	
	! Set radius (Stress Period 1)
	r(1)=1./12. ! Must NOT divide by zero, see below for vr
	
	do j=2,col
		r(j)=delr+r(j-1)
	end do

	! Set velocity as a funciton of radius, i.e., v(r)
	do j=1,col
		vr(j)=Qi/(2.*pi*r(j)*b*theta)
	end do
	
	! Compute
	! =============================================================================================
	! Stress period 1: forced-gradient injection
	! =============================================================================================
	
	nend=ti/deltS1
	time=0.
	vol=0.
	rmin=0.
	
	kount1=0
	kount2=nend/100
	eye=0

	! Zero out variables	
	do i=1,sim
		delM(i)=0.
		Disc(i)=0.
	end do
	
	! Explicit forward difference solution
	do n=1,nend
		do i=1,sim
			do j=2,col-1
				A=alpha(i)*((vr(j+1)+vr(j))/2.*(ColdS1(i,j+1)-ColdS1(i,j))/delr-(vr(j)+vr(j-1))/2.*&
				(ColdS1(i,j)-ColdS1(i,j-1))/delr)/delr-&
				(vr(j)+vr(j-1))/2.*(ColdS1(i,j)-ColdS1(i,j-1))/delr+k(i)*ColdS1(i,j)
				
				CnewS1(i,j)=A*deltS1/Ret(i)+ColdS1(i,j)
			end do
		end do
		
		! Mass discrepancy
		if(MB=="Yes" .or. MB=="yes") then
			! Reset masses
			do i=1,sim
				Mo(i)=0.
				Mt(i)=0.
				So(i)=0.
				Si(i)=0.
			end do
			! Old and new dissolved and sorbed (interior nodes)
			do i=1,sim
				do j=2,col-2
					Mo(i)=(((ColdS1(i,j+1)+ColdS1(i,j))/2.)*(r(j+1)**2.-r(j)**2.)*pi*b*theta*28.3168/1000.)*Ret(i)+&
					Mo(i)
					Mt(i)=(((CnewS1(i,j+1)+CnewS1(i,j))/2.)*(r(j+1)**2.-r(j)**2.)*pi*b*theta*28.3168/1000.)*Ret(i)+&
					Mt(i)
				end do
			end do
			do i=1,sim
				delM(i)=Mt(i)-Mo(i)+delM(i)
			end do	
			! External source and sink mass (inlet and outlet)
			do i=1,sim
				So(i)=((vr(1)+vr(2))/2.)*((ColdS1(i,1)+ColdS1(i,2)+CnewS1(i,1)+CnewS1(i,2))/4.)*deltS1*&
				pi*2.*b*theta*((r(1)+r(2))/2.)*28.3168/1000.
				Si(i)=((vr(col-1)+vr(col))/2.)*((ColdS1(i,col-1)+ColdS1(i,col)+CnewS1(i,col-1)+CnewS1(i,col))/4.)*deltS1*&
				pi*2.*b*theta*((r(col-1)+r(col))/2.)*28.3168/1000.
			end do
			! External sink mass (decay interior nodes)
			do i=1,sim
				do j=2,col-2
					Si(i)=((ColdS1(i,j+1)+ColdS1(i,j))/2.*(r(j+1)**2.-r(j)**2.)*pi*b*theta*28.3168/1000.)-&
					((ColdS1(i,j+1)+ColdS1(i,j))/2.*(r(j+1)**2.-r(j)**2.)*pi*b*theta*28.3168/1000.)*e**(k(i)*deltS1)+&					
					Si(i)
				end do
			end do
					
			! Cumulative mass discrepancy
			do i=1,sim
				Disc(i)=(((Mo(i)+So(i))-(Mt(i)+Si(i)))/(0.5*((Mo(i)+So(i))+(Mt(i)+Si(i)))))*100.+Disc(i)
			end do
		end if
		
		time=time+deltS1
		kount1=1+kount1
		
		! Output stress period 1 mass budget
		if(kount1==kount2) then
			kount1=0
			eye=1+eye
			DiscS1(eye,1)=time
			DiscS1(eye,2)=Disc(1)
			DiscS1(eye,3)=delM(1)
			DiscS1(eye,4)=Disc(2)
			DiscS1(eye,5)=delM(2)
			DiscS1(eye,6)=Disc(3)
			DiscS1(eye,7)=delM(3)
		end if
		
		! Reset old values (n) to new values (n+1)
		do j=2,col-1
			ColdS1(1,j)=CnewS1(1,j)
			ColdS1(2,j)=CnewS1(2,j)
			ColdS1(3,j)=CnewS1(3,j)
		end do
		
		! Keep track of time, volume, and radius, all increasing with time
		vol=Qi*time
		rmin=(vol/(pi*b*theta))**(1./2.)
		
	end do
	
	! Radial (1/12ft to 25ft) and Concs (sims 1, 2, & 3)
	do i=1,col
		Cr(i,1)=r(i)
		Cr(i,2)=CnewS1(1,i)
		Cr(i,3)=CnewS1(2,i)
		Cr(i,4)=CnewS1(3,i)
	end do

	! Horizontal
	do i=1,col+col-1
		Cx(i,1)=-25.0+(1./12.)*(i-1)
	end do
	
	! Right (+) side numerical
	do i=1,col
		Cx(i+col-1,2)=CnewS1(1,i)
		Cx(i+col-1,3)=CnewS1(2,i)
		Cx(i+col-1,4)=CnewS1(3,i)
	end do
	
	! Left (-) side numerical
	do i=2,col
		Cx(302-i,2)=CnewS1(1,i)
		Cx(302-i,3)=CnewS1(2,i)
		Cx(302-i,4)=CnewS1(3,i)
	end do
	
	! Analytical solution for 5, 6, and 7th column in C vs x simulations
	if (AM=="Yes" .or. AM=="yes") then
	do j=1,sim
		do i=1,col+col-1
			A=Cx(i,1)**2.-((Qi*ti)/(pi*b*theta))
			AA=((((((Qi*ti)/(pi*b*theta))**(3./2.))*16.*alpha(j)))/3.)**(1./2.)
				if(A/AA >= 0.) then
					AAA=1.-(1./(((1.+a1*(A/AA)+a2*(A/AA)**2.+a3*(A/AA)**3.+a4*(A/AA)**4.))**4.))
					AAA=1.-AAA
				else
					AAA=-1.*(1.-(1./(((1.+a1*(-A/AA)+a2*(-A/AA)**2.+a3*(-A/AA)**3.+a4*(-A/AA)**4.))**4.)))
					AAA=1.-AAA
				end if
			Cx(i,j+sim+1)=(Ci/2.)*AAA
		end do
	end do
	end if
	
	! Output Stres Period 1
	write(*,*) ""
	write(*,*) "Output: Injection Phase"
	write(*,*) ""
	HeaderS1(1)="#Distance (ft)"
	HeaderS1(2)="Sim. 1 (mg/L)"
	HeaderS1(3)="Sim. 2 (mg/L)"
	HeaderS1(4)="Sim. 3 (mg/L)"
	HeaderS1(5)="Ana. 1 (mg/L)"
	HeaderS1(6)="Ana. 2 (mg/L)"
	HeaderS1(7)="Ana. 3 (mg/L)"
	HeaderS1(:)=adjustr(HeaderS1(:))

	HeaderDiscS1(1)="#Time (hours)"
	HeaderDiscS1(2)="Disc. 1 (%)"
	HeaderDiscS1(3)="delM. 1 (grams)"
	HeaderDiscS1(4)="Disc. 2 (%)"
	HeaderDiscS1(5)="delM. 2 (grams)"
	HeaderDiscS1(6)="Disc. 3 (%)"
	HeaderDiscS1(7)="delM. 3 (grams)"
	HeaderDiscS1(:)=adjustr(HeaderDiscS1(:))

	! Command line output
	write(*,10) (HeaderS1(j),j=1,sim+1+sim)
	10 format(1x,7a20)
	write(*,11) ((Cx(i,j),j=1,sim+1+sim),i=1,col+col-1,5)
	11 format(1x,7f20.2)
	
	if(MB=="Yes" .or. MB=="yes") then
		write(*,*) ""
		write(*,*) "Mass Budget: Cumulative Discrepancy"
		write(*,*) ""
		write(*,12) Disc(1), Disc(2), Disc(3)
		12 format (1x,"Sim.1 = ",f10.3,"%",/,1x,"Sim.2 = ",f10.3,"%",/,1x,"Sim.3 = ",f10.3,"%")
	end if
	
	! ==========OPEN OUTPUT FILES==========
	open(unit=1000, file='swid_out_s1.conc') ! Open new or overwrite old output file for gnuplot (#) stress period 1
	open(unit=1001, file='swid_out_s2.conc') ! Open new or overwrite old output file for gnuplot (#) stress period 2
	open(unit=2000, file='swid_out_s1.disc') ! Open new or overwrite old output file for mass discrepancy stress period 1
	open(unit=2001, file='swid_out_s2.disc') ! Open new or overwrite old output file for mass discrepancy stress period 2
	
	write(1000,100) (HeaderS1(j),j=1,sim+1+sim)
	100 format(" #Output: Injection Phase",/1x,7a20)
	write(1000,101) ((Cx(i,j),j=1,sim+1+sim),i=1,col+col-1)
	101 format(1x,7f20.4)
	
	! User inputs singular
	write(1000,102) Qi/192.5, ti/0.04166667, b, theta, Ci, Ca
	102 format(" #Pumping rate (gpm)= ",f10.4,/" #Pumping time (hours)= ",f10.4,/" #Aquifer thickenss (ft)= ",f10.4,/&
	" #Aquifer porosity (-)= ",f10.4,/" #Solute concentration injection fluid (mg/L)= ",f10.4,/&
	" #Solute concentration aquifer fluid (mg/L)= ",f10.4)
	
	! User inputs of sims 1, 2, and 3	
	write(1000,103) (alpha(i),i=1,sim)
	103 format(" #Aquifer dispersivities for simulations 1, 2, & 3 (ft)= ",f10.4,",",f10.4,", &",f10.4)
	write(1000,104) (Ret(i),i=1,sim)
	104 format(" #Solute retentions for simulations 1, 2, & 3 (-)= ",f10.4,",",f10.4,", &",f10.4)
	write(1000,105) (k(i),i=1,sim)
	105 format(" #Solute decay rates for simulations 1, 2, & 3 (1/day)= ",f10.4,",",f10.4,", &",f10.4)

	! Mass discrepancy
	if(MB=="Yes" .or. MB=="yes") then
		do i=1,100
			DiscS1(i,1)=DiscS1(i,1)/0.04166667 ! Convert time in days to time in hours
		end do
		write(2000,107) (HeaderDiscS1(j),j=1,sim+1+sim)
		107 format(" #Mass discrepancy: Injection Phase",/1x,7a20)
		write(2000,108) ((DiscS1(i,j),j=1,sim+1+sim),i=1,100)
		108 format(1x,7f20.6)
	end if

	! Compute
	! =============================================================================================
	! Stress period 2: natural-gradient drift
	! =============================================================================================
	
	! Set boundary conditions at inlet (left side) and outlet (right side)
	do i=1,sim
		ColdS2(i,1)=Ca ! Inlet
		CnewS2(i,1)=Ca ! Inlet
	end do
	
	do i=1,sim
		ColdS2(i,col+col-1)=Ca ! Outlet
		CnewS2(i,col+col-1)=Ca ! Outlet
	end do
	
	! Set initial condition within boundaries
	do j=2,col+col-1-1
		ColdS2(1,j)=Cx(j,2)
		ColdS2(2,j)=Cx(j,3)
		ColdS2(3,j)=Cx(j,4)
	end do
	
	nend=30000
	kount1=0
	kount2=100
	time=0.
	i=0 ! This "i" can be problem...use "eye" as a work around when needed
	
	kount3=0
	kount4=nend/100
	aye=0
	
	! Zero out variables
	do eye=1,sim
		delM(eye)=0.
		Disc(eye)=0.
	end do
	
	write(*,*) ""
	write(*,*) "Output: Drift Phase"
	write(*,*) ""
	HeaderS2(1)="#Time (days)"
	HeaderS2(2)="Sim.1 (mg/L)"
	HeaderS2(3)="Sim.2 (mg/L)"
	HeaderS2(4)="Sim.3 (mg/L)"
	HeaderS2(:)=adjustr(HeaderS2(:))
	
	HeaderDiscS2(1)="#Time (days)"
	HeaderDiscS2(2)="Disc. 1 (%)"
	HeaderDiscS2(3)="delM. 1 (grams)"
	HeaderDiscS2(4)="Disc. 2 (%)"
	HeaderDiscS2(5)="delM. 2 (grams)"
	HeaderDiscS2(6)="Disc. 3 (%)"
	HeaderDiscS2(7)="delM. 3 (grams)"
	HeaderDiscS2(:)=adjustr(HeaderDiscS2(:))
	
	! Explicit forward difference solution
	do n=1,nend
		do eye=1,sim
			do j=2,col+col-1-1
				A=vx(eye)*alpha(eye)*((ColdS2(eye,j+1)-2.*ColdS2(eye,j)+ColdS2(eye,j-1))/(delx**2.))-&
				vx(eye)*((ColdS2(eye,j)-ColdS2(eye,j-1))/(delx))+k(eye)*ColdS2(eye,j)
				CnewS2(eye,j)=A*deltS2/Ret(eye)+ColdS2(eye,j)
			end do
		end do
		
		! Mass discrepancy
		if(MB=="Yes" .or. MB=="yes") then
			! Reset masses
			do eye=1,sim
				Mo(eye)=0.
				Mt(eye)=0.
				So(eye)=0.
				Si(eye)=0.
			end do
			! Old and new dissolved and sorbed mass (interior nodes)
			do eye=1,sim
				do j=2,col+col-1-1
					Mo(eye)=((ColdS2(eye,j+1)+ColdS2(eye,j))/2.*delx*delx*b*theta*28.3168/1000.)*Ret(eye)+&
					Mo(eye)
					Mt(eye)=((CnewS2(eye,j+1)+CnewS2(eye,j))/2.*delx*delx*b*theta*28.3168/1000.)*Ret(eye)+&
					Mt(eye)
				end do
			end do
			do eye=1,sim
				delM(eye)=Mt(eye)-Mo(eye)+delM(eye)
			end do
			! External source and sink mass (inlet and outlet)
			do eye=1,sim
				So(eye)=vx(eye)*((ColdS2(eye,1)+ColdS2(eye,2)+CnewS2(eye,1)+CnewS2(eye,2))/4.)*deltS2*delx*b*theta*28.3168/1000.			
				Si(eye)=vx(eye)*((ColdS2(eye,col+col-1)+ColdS2(eye,col+col-2)+CnewS2(eye,col+col-1)+CnewS2(eye,col+col-2))/4.)&
				*deltS2*delx*b*theta*28.3168/1000.
			end do
			! External sink mass (decay interior nodes)
			do eye=1,sim
				do j=2,col+col-1-1
					Si(eye)=((ColdS2(eye,j+1)+ColdS2(eye,j))/2.*delx*delx*b*theta*28.3168/1000.)-&
					((ColdS2(eye,j+1)+ColdS2(eye,j))/2.*delx*delx*b*theta*28.3168/1000.)*e**(k(eye)*deltS2)+&					
					Si(eye)
				end do
			end do
			! Cumulative mass discrepancy
			do eye=1,sim
				Disc(eye)=(((Mo(eye)+So(eye))-(Mt(eye)+Si(eye)))/(0.5*((Mo(eye)+So(eye))+(Mt(eye)+Si(eye)))))*100.
			end do
		end if
		
		! Reset old values (n) to new values (n+1)
		do j=2,col+col-1-1
			ColdS2(1,j)=CnewS2(1,j)
			ColdS2(2,j)=CnewS2(2,j)
			ColdS2(3,j)=CnewS2(3,j)
		end do
		
		time=time+deltS2
		kount1=1+kount1
		kount3=1+kount3	
		
	! Output stress period 2
		if(kount1==kount2) then
			kount1=0
			i=1+i
			Ct(i,1)=time
			Ct(i,2)=CnewS2(1,injw)
			Ct(i,3)=CnewS2(2,injw)
			Ct(i,4)=CnewS2(3,injw)
		end if
		
		if (kount3==kount4) then
			kount3=0
			aye=1+aye
			DiscS2(aye,1)=time
			DiscS2(aye,2)=Disc(1)
			DiscS2(aye,3)=delM(1)
			DiscS2(aye,4)=Disc(2)
			DiscS2(aye,5)=delM(2)
			DiscS2(aye,6)=Disc(3)
			DiscS2(aye,7)=delM(3)
		end if
	end do

	! Command line output
	write(*,20) (HeaderS2(j),j=1,sim+1)
	20 format(1x,4a20)
	write(*,21) ((Ct(i,j),j=1,sim+1),i=1,col-1,5)
	21 format(1x,4f20.2)
	
	if(MB=="Yes" .or. MB=="yes") then
		write(*,*) ""
		write(*,*) "Mass Budget: Cumulative Discrepancy"
		write(*,*) ""
		write(*,22) Disc(1), Disc(2), Disc(3)
		22 format (1x,"Sim.1 = ",f10.3,"%",/,1x,"Sim.2 = ",f10.3,"%",/,1x,"Sim.3 = ",f10.3,"%")
	end if

	write(1001,200) (HeaderS2(j),j=1,sim+1)
	200 format(" #Output: Drift Phase",/1x,4a20)
	write(1001,201) ((Ct(i,j),j=1,sim+1),i=1,col-1)
	201 format(1x,4f20.4)
	
	! User inputs of sims 1, 2, and 3	
	write(1001,202) (vx(i),i=1,sim)
	202 format(" #Groundwater velocities for simulations 1, 2, & 3 (ft/day)= ",f10.4,",",f10.4,", &",f10.4)
	write(1001,203) (alpha(i),i=1,sim)
	203 format(" #Aquifer dispersivities for simulations 1, 2, & 3 (ft)= ",f10.4,",",f10.4,", &",f10.4)
	write(1001,204) (Ret(i),i=1,sim)
	204 format(" #Solute retentions for simulations 1, 2, & 3 (-)= ",f10.4,",",f10.4,", &",f10.4)
	write(1001,205) (k(i),i=1,sim)
	205 format(" #Solute decay rates for simulations 1, 2, & 3 (1/day)= ",f10.4,",",f10.4,", &",f10.4)
	
	! Mass discrepancy
	if(MB=="Yes" .or. MB=="yes") then
		!do i=1,100
		!	DiscS2(i,1)=DiscS2(i,1)/0.04166667 ! Convert time in days to time in hours
		!end do
		write(2001,206) (HeaderDiscS2(j),j=1,sim+1+sim)
		206 format(" #Mass discrepancy: Drift Phase",/1x,7a20)
		write(2001,207) ((DiscS2(i,j),j=1,sim+1+sim),i=1,100)
		207 format(1x,7f20.6)
	end if
	
	! ==========CLOSE OUTPUT FILES==========
	close(1000)
	close(1001)
	close(2000)
	close(2001)

end program swid_com_23