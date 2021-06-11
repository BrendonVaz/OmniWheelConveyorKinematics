
import cv2 					as cv
import numpy 				as np
import random 				as rd
import math 				as mt
import matplotlib
import matplotlib.pyplot 	as pl
import matplotlib.patches 	as pa
import numpy.linalg			as la
import math 				as mt
import time

def main(args):
	
	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#~ constants
	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	cellangl 	= [  0,180,  0,180,  0,180,180,  0,180,  0,180,  0,  0,180,  0,180,  0,180,180,  0,180,  0,180,  0,  0,180,  0,180,  0,180,180,  0,180,  0,180,  0]; 	#angle of each cell

	umx 		= 15.;																																					#maximum allowed wheel speed
	umn 		= 2.;																																					#maximum allowed wheel speed
	tau 		= .2;																																					#time step
	xyt 		= np.array([15.0,15.0, 0*mt.pi/180])																													#target		final
	est 		= np.array([15.0,15.0, 0*mt.pi/180])																													#estimate	initial
	pl.ion();

	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#~ set up & layout
	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	k=0; motopos = {}; senspos = {}; motovel = {};																														#calculate installation geometry

	for i in range(1,7):
		for j in range(1,7):
			motopos[(i,j,1)] = ((i-1)*5.0+2.5 + (-1.25)*mt.cos(cellangl[k]*mt.pi/180) - (-1.25)*mt.sin(cellangl[k]*mt.pi/180),(j-1)*5.0+2.5 + (-1.25)*mt.sin(cellangl[k]*mt.pi/180) + (-1.25)*mt.cos(cellangl[k]*mt.pi/180) , (0)*mt.cos(cellangl[k]*mt.pi/180) - (1)*mt.sin(cellangl[k]*mt.pi/180), (0)*mt.sin(cellangl[k]*mt.pi/180) + (1)*mt.cos(cellangl[k]*mt.pi/180));
			motopos[(i,j,2)] = ((i-1)*5.0+2.5 + ( 1.25)*mt.cos(cellangl[k]*mt.pi/180) - ( 1.25)*mt.sin(cellangl[k]*mt.pi/180),(j-1)*5.0+2.5 + ( 1.25)*mt.sin(cellangl[k]*mt.pi/180) + ( 1.25)*mt.cos(cellangl[k]*mt.pi/180) , (1)*mt.cos(cellangl[k]*mt.pi/180) - (0)*mt.sin(cellangl[k]*mt.pi/180), (1)*mt.sin(cellangl[k]*mt.pi/180) + (0)*mt.cos(cellangl[k]*mt.pi/180));
			senspos[(i,j,1)] = ((i-1)*5.0+2.5 + (-1.25)*mt.cos(cellangl[k]*mt.pi/180) - ( 1.25)*mt.sin(cellangl[k]*mt.pi/180),(j-1)*5.0+2.5 + (-1.25)*mt.sin(cellangl[k]*mt.pi/180) + ( 1.25)*mt.cos(cellangl[k]*mt.pi/180));
			senspos[(i,j,2)] = ((i-1)*5.0+2.5 + ( 1.25)*mt.cos(cellangl[k]*mt.pi/180) - (-1.25)*mt.sin(cellangl[k]*mt.pi/180),(j-1)*5.0+2.5 + ( 1.25)*mt.sin(cellangl[k]*mt.pi/180) + (-1.25)*mt.cos(cellangl[k]*mt.pi/180));
			k+=1;

	for i in range(1,7):
		for j in range(1,7):
			for k in range(1,3):
				if abs(motopos[(i,j,k)][2]) < 0.0001: a = list(motopos[(i,j,k)]);a[2] = 0.0;motopos[(i,j,k)] = tuple(a);
				if abs(motopos[(i,j,k)][3]) < 0.0001: a = list(motopos[(i,j,k)]);a[3] = 0.0;motopos[(i,j,k)] = tuple(a);
				print motopos[(i,j,k)]

	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#~ define state kinematics
	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	#~ let v_dot = [v1, v2, v3... v72] = velocity from wheel
	#~ let x_dot = [xdot, ydot, tdot]  = change in state x,y,t where (x,y) = position and (t) = orientation of material...

	#~ from geometry/kinematics, we know that for each wheel i we can formulate row matrix [A]
	
	#~ v_dot = ( [mu] [mv] [mv*(mx-cx) - mu*(my-cy)] ) dot (xdot, ydot, tdot) 
	#~ v_dot = ( [A] ) dot (xdot, ydot, tdot)...  where row [A]_i = ( [mu] [mv] [mv*(mx-x) - mu*(my-y)] ) for each wheel i
	
	#~ taking the pseudoinverse of A we arrive at state kinematics given input v_dot

	#~ x_dot = [B] dot (v_dot, 1 to n), where [B] = inv([A]'[A])*[A]... [A]'  = Transpose of A 

	#~ notice dimensions - 3xn * nx3 * 3xn = 3xn x nx1 = 3x1 vector of state quantities

	#~ mu,mv	= velocity vector of motor i:			known from geometry
	#~ mx,my	= position vector of motor i:			known from geometry
	#~ cx,cy	= position vector of material:			known from current state of X
	
	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#~ define error vector
	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	#~ xyt 		= [x_d, y_d, t_d]	= state desired  	= Next point in path to follow 
	#~ est	 	= [x_c, y_c, t_c]	= state estimate 	= solution from kalman filter with input [mbr algorithm] + [propogated state]
	#~ err 		= [x_e, y_e, t_e]	= state error		= designed state - designed state

	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#~ define control law u = k_p*u(err) + k_i*int(u(err)) + kd*der(u(err)), v_dot = -kp*u(x); where kp=1 for proportional feedback control with gain of 1 unit
	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	#~ vdot = ( [mu] [mv] [mv*(mx-cx) - mu*(my-cy)] ) dot (x_e, y_e, t_e)

	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	xyt = np.array([15.0,15.0, 90*mt.pi/180])						#state target
	est = np.array([15.0,15.0, 90*mt.pi/180])						#state estimate
	err = np.array([ 1.0, 1.0,  1*mt.pi/180])						#state error
	
	if (-90*mt.pi/180-est[2]*mt.pi/180) < (90*mt.pi/180-est[2]*mt.pi/180): xyt[2] = -90*mt.pi/180

	#~ for time_step in range(0,25):

	time_step = 0;state = 0;


	try:
		while True:

			time_step+=1;													#state error
			err = xyt-est													#state error
			dxi = err*0														#state update
			
			if err[0] < 0.05 and err[1] < 0.05 and err[2] < 5*mt.pi/180:
				if state==1: xyt = np.array([30.0,15.0, xyt[2]]);state=2; err = np.array([ 1.0, 1.0,  1*mt.pi/180])
				else:		 xyt = np.array([15.0,15.0, xyt[2]]);state=1; est = np.array([rd.uniform(10,20),rd.uniform(10,20),rd.uniform(-45,45)*mt.pi/180]);pl.figure(1);pl.cla();pl.figure(2);pl.cla();

			if (-90*mt.pi/180-est[2]*mt.pi/180) < (90*mt.pi/180-est[2]*mt.pi/180): xyt[2] = -90*mt.pi/180		#we want vertical, dont care about facing left/right, -90/90

			a = np.empty((72,  3))*0.;										# inverse kinematics matrix
			b = np.empty(( 3, 72))*0.;										# inverse kinematics matrix transpose
			c = np.empty(( 3,  3))*0.;										# A.T * A
			d = np.empty(( 3,  3))*0.;										# inv(A.T * A)
			e = np.empty(( 3, 72))*0.;										# forward kinematics matrix result
			f = np.empty(( 3,  3))*0.;										# verify f = E*A = Identity
			u = np.empty((72,  1))*0.;										# control input K*U(x)

			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			#~ inverse kinematics & control input
			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

			m = 0;
			for j in range(1,7): #motovel for display only... 
				for i in range(1,7):
					for k in range(1,3):
						motovel[i,j,k] = motopos[i,j,k][2]*err[0]  +  motopos[i,j,k][3]*err[1]  +  (-(motopos[i,j,k][1]-est[1])*err[2]*motopos[i,j,k][2] + (motopos[i,j,k][0]-est[0])*err[2]*motopos[i,j,k][3]); 
			m=0;
			for i in range(1,7):
				for j in range(1,7):
					for k in range(1,3):
						a[m][0] = motopos[i,j,k][2]
						a[m][1] = motopos[i,j,k][3]				
						a[m][2] = motopos[i,j,k][3]*(motopos[i,j,k][0]-est[0]) - motopos[i,j,k][2]*(motopos[i,j,k][1]-est[1])
							
						u[m][0] = motopos[i,j,k][2]*err[0]  +  motopos[i,j,k][3]*err[1]  +  (-(motopos[i,j,k][1]-est[1])*err[2]*motopos[i,j,k][2] + (motopos[i,j,k][0]-est[0])*err[2]*motopos[i,j,k][3]); 
						m+=1;

			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			#~ control input saturation betweeu u_max and u_min if u>0
			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

			if len(u[u!=0]) > 0:

				smx = min(umx,max(abs(u)))
				smn = max(umn,min(abs(u)))
				
				u[u>0] =  ((smx-umn)*(abs(u[u>0])-min(abs(u[u>0])))/(max(abs(u[u>0]))-min(abs(u[u>0]))) + umn)
				u[u<0] = -((smx-umn)*(abs(u[u<0])-min(abs(u[u<0])))/(max(abs(u[u<0]))-min(abs(u[u<0]))) + umn)

			a[np.abs(a) < 1e-10] = 0.; 
			
			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			#~ state propogation
			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			
			b = np.transpose(a)		;	b[np.abs(b) < 1e-10] = 0.
			c = np.matmul(b,a)		; 	c[np.abs(c) < 1e-10] = 0.
			d = la.inv(c)			;	d[np.abs(d) < 1e-10] = 0.
			e = np.matmul(d,b)		;	e[np.abs(e) < 1e-10] = 0.
			f = np.matmul(e,a)		;	f[np.abs(f) < 1e-10] = 0.
			
			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			#~ state & error update
			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

			print "------------------------------------------------\n", "State Target:\n",			"------------------------------------------------\n", xyt
			print "------------------------------------------------\n", "State Error:\n",			"------------------------------------------------\n", err
			print "------------------------------------------------\n", "State Estimate:\n",		"------------------------------------------------\n", est

			dxi	= np.matmul(e,u)*tau;
			est = np.reshape(est+np.transpose(dxi),(3,))

			print "------------------------------------------------\n", "Update:\n",				"------------------------------------------------\n", np.transpose(dxi)
			print "------------------------------------------------\n", "New Estimate:\n",			"------------------------------------------------\n", est
			print "------------------------------------------------\n", "New Error:\n",				"------------------------------------------------\n", xyt-est
			print "------------------------------------------------\n", "Elapsed Time:\n",			"------------------------------------------------\n", (time_step+1)*tau

			print "------------------------------------------------\n", "Inverse Kinematics:\n",	"------------------------------------------------\n", a
			print "------------------------------------------------\n", "Forward Kinematics:\n",	"------------------------------------------------\n", e
			print "------------------------------------------------\n", "Control Input:\n",			"------------------------------------------------\n", u

			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			#~ Dyanmic visualization
			#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

			pl.figure(1)

			pl.plot(time_step,err[0],'.-r');
			pl.plot(time_step,err[1],'.-g');
			pl.plot(time_step,err[2],'.-b');
			pl.pause(0.01)


			pl.figure(2)

			mx=[];my=[];mu=[];mv=[];ru=[];rv=[];sx=[];sy=[];m = 0;

			cor = (12,9)

			crx = []; cry = [];	#core vertices current
			cfx = []; cfy = []; #core vertices target

			for i in range(1,7):
				for j in range(1,7):
					for k in range(1,3):
						mx.append(motopos[i,j,k][0]);my.append(motopos[i,j,k][1]);
						mu.append(motopos[i,j,k][2]);mv.append(motopos[i,j,k][3]);
						sx.append(senspos[i,j,k][0]);sy.append(senspos[i,j,k][1]);
						ru.append(motovel[i,j,k] * motopos[i,j,k][2]);
						rv.append(motovel[i,j,k] * motopos[i,j,k][3]);

			cfx.append(xyt[0] + ( cor[0])*(mt.cos(xyt[2])) - ( cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + ( cor[0])*(mt.sin(xyt[2])) + ( cor[1])*(mt.cos(xyt[2])));
			cfx.append(xyt[0] + ( cor[0])*(mt.cos(xyt[2])) - (-cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + ( cor[0])*(mt.sin(xyt[2])) + (-cor[1])*(mt.cos(xyt[2])));
			cfx.append(xyt[0] + (-cor[0])*(mt.cos(xyt[2])) - (-cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + (-cor[0])*(mt.sin(xyt[2])) + (-cor[1])*(mt.cos(xyt[2])));
			cfx.append(xyt[0] + (-cor[0])*(mt.cos(xyt[2])) - ( cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + (-cor[0])*(mt.sin(xyt[2])) + ( cor[1])*(mt.cos(xyt[2])));
			cfx.append(xyt[0] + ( cor[0])*(mt.cos(xyt[2])) - ( cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + ( cor[0])*(mt.sin(xyt[2])) + ( cor[1])*(mt.cos(xyt[2])));

			crx.append(est[0] + ( cor[0])*(mt.cos(est[2])) - ( cor[1])*(mt.sin(est[2]))); cry.append(est[1] + ( cor[0])*(mt.sin(est[2])) + ( cor[1])*(mt.cos(est[2])));
			crx.append(est[0] + ( cor[0])*(mt.cos(est[2])) - (-cor[1])*(mt.sin(est[2]))); cry.append(est[1] + ( cor[0])*(mt.sin(est[2])) + (-cor[1])*(mt.cos(est[2])));
			crx.append(est[0] + (-cor[0])*(mt.cos(est[2])) - (-cor[1])*(mt.sin(est[2]))); cry.append(est[1] + (-cor[0])*(mt.sin(est[2])) + (-cor[1])*(mt.cos(est[2])));
			crx.append(est[0] + (-cor[0])*(mt.cos(est[2])) - ( cor[1])*(mt.sin(est[2]))); cry.append(est[1] + (-cor[0])*(mt.sin(est[2])) + ( cor[1])*(mt.cos(est[2])));
			crx.append(est[0] + ( cor[0])*(mt.cos(est[2])) - ( cor[1])*(mt.sin(est[2]))); cry.append(est[1] + ( cor[0])*(mt.sin(est[2])) + ( cor[1])*(mt.cos(est[2])));

			pl.cla(); #comment to show path
			
			pl.plot(mx,my,'o', markersize=4,color='g')											#motor  position marker
			pl.quiver(mx,my,ru,rv,color='b',  headwidth=2,headlength=3,width=0.004) 			#motor velocity ctrl vector
			pl.plot(crx,cry,'.-' , color='gold')
			pl.plot(cfx,cfy,'.--', color='black')
			
			pl.gca().set_aspect('equal',adjustable='box');pl.xlim(-5,35); pl.ylim(-5,35);pl.pause(0.01)

	except KeyboardInterrupt: sys.exit();

	pl.ioff();pl.show();
	
	return












	#~ for i in range(20):
		#~ a+=i%1;
		#~ b+=i%2;
		#~ c+=i%3;
		#~ d+=i%4;
		#~ e+=i%5;
		#~ pl.cla();pl.gca().set_aspect('equal',adjustable='box');pl.xlim(-5,35); pl.ylim(-5,35);pl.ion()				#create plot...
		#~ pl.figure(1);pl.plot(a,x,'or')
		#~ pl.figure(1);pl.plot(b,x,'og')
		#~ pl.figure(1);pl.plot(c,x,'ob')
		#~ pl.figure(1);pl.plot(d,x,'oc')
		#~ pl.figure(1);pl.plot(e,x,'om')

		#~ pl.pause(0.1)
		
	#~ pl.ioff();pl.show()

	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	#~ static visualization
	#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	#~ for j in range(1,7):
		#~ for i in range(1,7):
			#~ for k in range(1,3):
				#~ motovel[i,j,k] = motopos[i,j,k][2]*err[0]  +  motopos[i,j,k][3]*err[1]  +  (-(motopos[i,j,k][1]-est[1])*err[2]*motopos[i,j,k][2] + (motopos[i,j,k][0]-est[0])*err[2]*motopos[i,j,k][3]); 

	#~ mx=[];my=[];mu=[];mv=[];ru=[];rv=[];sx=[];sy=[];

	#~ for j in range(1,7):
		#~ for i in range(1,7):
			#~ for k in range(1,3):
				#~ mx.append(motopos[i,j,k][0]);my.append(motopos[i,j,k][1]);
				#~ mu.append(motopos[i,j,k][2]);mv.append(motopos[i,j,k][3]);
				#~ sx.append(senspos[i,j,k][0]);sy.append(senspos[i,j,k][1]);
				#~ ru.append(motovel[i,j,k] * motopos[i,j,k][2]);
				#~ rv.append(motovel[i,j,k] * motopos[i,j,k][3]);

 	#~ cor = (12,9)

	#~ crx = []; cry = [];	#core vertices current
	#~ cfx = []; cfy = []; #core vertices target

	#~ cfx.append(xyt[0] + ( cor[0])*(mt.cos(xyt[2])) - ( cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + ( cor[0])*(mt.sin(xyt[2])) + ( cor[1])*(mt.cos(xyt[2])));
	#~ cfx.append(xyt[0] + ( cor[0])*(mt.cos(xyt[2])) - (-cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + ( cor[0])*(mt.sin(xyt[2])) + (-cor[1])*(mt.cos(xyt[2])));
	#~ cfx.append(xyt[0] + (-cor[0])*(mt.cos(xyt[2])) - (-cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + (-cor[0])*(mt.sin(xyt[2])) + (-cor[1])*(mt.cos(xyt[2])));
	#~ cfx.append(xyt[0] + (-cor[0])*(mt.cos(xyt[2])) - ( cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + (-cor[0])*(mt.sin(xyt[2])) + ( cor[1])*(mt.cos(xyt[2])));
	#~ cfx.append(xyt[0] + ( cor[0])*(mt.cos(xyt[2])) - ( cor[1])*(mt.sin(xyt[2]))); cfy.append(xyt[1] + ( cor[0])*(mt.sin(xyt[2])) + ( cor[1])*(mt.cos(xyt[2])));

	#~ crx.append(est[0] + ( cor[0])*(mt.cos(est[2])) - ( cor[1])*(mt.sin(est[2]))); cry.append(est[1] + ( cor[0])*(mt.sin(est[2])) + ( cor[1])*(mt.cos(est[2])));
	#~ crx.append(est[0] + ( cor[0])*(mt.cos(est[2])) - (-cor[1])*(mt.sin(est[2]))); cry.append(est[1] + ( cor[0])*(mt.sin(est[2])) + (-cor[1])*(mt.cos(est[2])));
	#~ crx.append(est[0] + (-cor[0])*(mt.cos(est[2])) - (-cor[1])*(mt.sin(est[2]))); cry.append(est[1] + (-cor[0])*(mt.sin(est[2])) + (-cor[1])*(mt.cos(est[2])));
	#~ crx.append(est[0] + (-cor[0])*(mt.cos(est[2])) - ( cor[1])*(mt.sin(est[2]))); cry.append(est[1] + (-cor[0])*(mt.sin(est[2])) + ( cor[1])*(mt.cos(est[2])));
	#~ crx.append(est[0] + ( cor[0])*(mt.cos(est[2])) - ( cor[1])*(mt.sin(est[2]))); cry.append(est[1] + ( cor[0])*(mt.sin(est[2])) + ( cor[1])*(mt.cos(est[2])));
	
	#~ pl.gca().set_aspect('equal',adjustable='box');pl.xlim(-5,35); pl.ylim(-5,35);						#plot formatting... 
	
	#~ pl.figure(1);	pl.plot(sx,sy,'o', markersize=4,color='r')											#sensor position marker
	#~ pl.figure(1);	pl.plot(mx,my,'o', markersize=4,color='g')											#motor  position marker

	#~ pl.figure(1);	pl.quiver(mx,my,ru,rv,color='b',  headwidth=2,headlength=3,width=0.004) 			#motor velocity ctrl vector

	#~ pl.figure(1);	pl.plot(crx,cry,'.-' , color='gold')												#core position estimate
	#~ pl.figure(1);	pl.plot(cfx,cfy,'.--', color='gold')												#core position target

	#~ pl.show()


	#~ time.sleep(1)
	#~ pl.cla()
	#~ time.sleep(1)
	#~ pl.close()
	
if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
