#include<iostream>
#include<stdio.h>
#include<cmath>
#define Num_Max 50
#define tol 0.1
//====================================================================================================================
//====================================================================================================================
int IMAX, JMAX;
double Lx, Ly;
double vin;
double dt1;

double rho_g;				// Density of gas phase
double mu_g;				// Viscosity of the gas phase
double c_mu;				// Constant for determination of eddy viscosity
double c_1, c_2;				
double sigma_k, sigma_e;
	
double x[Num_Max], y[Num_Max];
double dx[Num_Max], dy[Num_Max];

int State[Num_Max][Num_Max];
double u_left, u_right;
double u_top, u_bottom;
double v_left, v_right;
double v_top, v_bottom;

double Ug[Num_Max][Num_Max];
double f[Num_Max][Num_Max];
int Conv_U;
double time_Ugb,time_Uga, cpu_time_U_gas;

double Vg[Num_Max][Num_Max];
double g[Num_Max][Num_Max];
int Conv_V;
double time_Vgb,time_Vga, cpu_time_V_gas;

double P[Num_Max][Num_Max];
int  Conv_P;

double KE[Num_Max][Num_Max];
double ke1[Num_Max][Num_Max];
int  Conv_KE;

double DE[Num_Max][Num_Max];
double de1[Num_Max][Num_Max];
int Conv_DE;
double solution_time;
//====================================================================================================================
// User-defined Header files
double sqr(double);
double mut(int, int);
double mue(int, int);
void input_data(void);
void grid_generation(void);
void domain_conditions(void);
void initialize(void);
void enforce_physics(void);
void u_update(void);
void v_update(void);
void P_update(void);
void correction(void);
void KE_update(void);
void DE_update(void);
void Steady_State_report(void);
//====================================================================================================================
using namespace std;
int main ()
{
	input_data();
	grid_generation();
	domain_conditions();
	initialize();
	enforce_physics();
	
	solution_time = 0.00;
	for (int nstep = 0; nstep < 10000; nstep++)
	{
		solution_time = dt1*nstep;
		if( (nstep%100) == 0 ){
			printf("Solution time  :: %lf\n", solution_time);
			//Transient_report(nstep);
		}
		u_update();
		v_update();
		P_update();
		correction();
		enforce_physics();
		KE_update();
		DE_update();
		enforce_physics();
		Steady_State_report();
	}	
	return 0;
}
//----------------------------------------------------------------------
void input_data()
{
	// Domain Information
	Lx = 5.0;
	Ly = 5.0;
	IMAX = 10;
	JMAX = 10;
	dt1 = 1e-6;
	// Gas Properties 
	rho_g	= 1.178;
	mu_g	= 1.983e-5;
	// Turbulence Model
	c_mu	= 0.09;
	c_1		= 1.44;
	c_2		= 1.92;
	sigma_k	= 1.00;
	sigma_e = 1.22;
	// Domain Wall info
	u_left	= 0.00;
	u_right	= 0.00;
	v_top	= 0.00;
	v_bottom= 0.00;
}
//----------------------------------------------------------------------
void grid_generation()
{
	int i, j, k;	
	printf("Grid Resolution :: IMAX=%d x JMAX=%d \n",IMAX, JMAX);
	printf("Grid Type :: UNIFORM CARTESIAN GRID\n");
	printf("---------------------------------------------------\n\n");
	//------------x-direction
	x[0] = - 0.50 * Lx/(IMAX);
	for (int i = 1; i <= IMAX; i++)
	{
		dx[i] = Lx/(IMAX);
		if(i == 1)
			x[i] = 0.50 * dx[i];
		else
			x[i] = x[i-1] + 0.50 * (dx[i] + dx[i-1]);
	} 
	x[IMAX+1] = Lx + 0.50 * Lx/(IMAX);	
	//------------y-direction
	y[0]= - 0.50 * Ly/(JMAX);
	for (int j = 1; j <= JMAX; j++)
	{
		dy[j] = Ly/(JMAX);
		if(j == 1) 
			y[j] = 0.50 * dy[j];
		else 
			y[j] = y[j-1] + 0.50 * (dy[j] + dy[j-1]);
	}
	y[JMAX+1] = Ly + 0.50 * Ly/(JMAX);  
}
//----------------------------------------------------------------------
void domain_conditions()
{
	//------------Defining Initial State
	for(int i = 0; i<=IMAX+1;i++)
	{
		for(int j = 0; j<=JMAX+1; j++)
		{
			State[i][j] = 0;
		}
	}
	
	for(int j = 0; j<=JMAX+1; j++)
	{
		if (State[0][j] == 0)		State[0][j]		= 1;					// Left Surface
	}
	
	//------------Defining Boundary Walls
	for(int j = 0; j<=JMAX+1; j++)
	{
		if (State[0][j] == 0)		State[0][j]		= -1;					// Left Surface
		if (State[IMAX+1][j] == 0)	State[IMAX+1][j] = -1;					// Right Surface
	}
	
	for(int i = 0; i<=IMAX+1;i++)
	{			
		if(State[i][0] == 0)		State[i][0]		= -1;					// Bottom Surface
		if(State[i][JMAX+1] == 0)	State[i][JMAX+1] = -1;					// Top Surface
	}
}
//----------------------------------------------------------------------
void initialize()
{
	vin = 10.0;
	double lm = 0.005;
	for(int i = 0; i<=IMAX+1; i++)
	{
		for(int j = 0; j<=JMAX+1; j++)
		{		
			Ug[i][j]	= 0.00;
			Vg[i][j]	= 0.00;
			P[i][j]		= 0.00;
			KE[i][j]	= 0.009*sqr(vin);
			DE[i][j]	= pow(0.009*sqr(vin),1.5)/lm;
		
			f[i][j]		= Ug[i][j];
			g[i][j]		= Vg[i][j];
			ke1[i][j]	= KE[i][j];
			de1[i][j]	= DE[i][j];
		}
	}
}
//----------------------------------------------------------------------
void enforce_physics()
{
	for(int j = 0; j<=JMAX+1; j++)
	{
		Ug[0][j]		= u_left;
		Vg[0][j]		= v_left;
		P[0][j]			= 0.00;
		
		Ug[IMAX+1][j]	= u_right;
		Vg[IMAX+1][j]	= v_right;
		P[IMAX+1][j]	= 0.00;
	}
	
	for(int i = 0; i<=IMAX+1;i++)
	{			
		Ug[i][0]		= u_bottom;
		Vg[i][0]		= v_bottom;
		P[i][0]			= 0.00;
		
		Ug[i][JMAX+1]	= u_top;
		Vg[i][JMAX+1]	= v_top;
		P[i][JMAX+1]	= 0.00;
	}	

	for(int i = 1; i <= IMAX; i++)
	{
		for(int j = 1; j <= JMAX; j++)
		{
			if (State[i][j] == -1)				// INTERIOR OBSTACLES REGIONS
			{				
				Ug[i][j]= 0.00;
				Vg[i][j]= 0.00;				
				f[i][j]	= Ug[i][j];
				g[i][j]	= Vg[i][j];
			}
			
			if (State[i][j] == 1)			// INLET REGIONS
			{				
				Ug[i][j]= vin;
				Vg[i][j]= 0.00;
				P[i][j]	= 0.00;
				
				f[i][j]	= Ug[i][j];
				g[i][j]	= Vg[i][j];
			}
			else if (State[i][j] == 2)		// OUTLET REGIONS
			{				
				P[i][j] = 0.00;
			}				
		}
	}
}
//----------------------------------------------------------------------
void u_update()
{
	//------Local variables declaration
	double bound = 100.0;
	double Fe_x[IMAX+1][JMAX+1], Fw_x[IMAX+1][JMAX+1];
	double Fn_x[IMAX+1][JMAX+1], Fs_x[IMAX+1][JMAX+1];
	double De_x[IMAX+1][JMAX+1], Dw_x[IMAX+1][JMAX+1];
	double Dn_x[IMAX+1][JMAX+1], Ds_x[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1];
	double Peclet_y[IMAX+1][JMAX+1];
	double Source_x[IMAX+1][JMAX+1];
	double Previous[IMAX+1][JMAX+1];
	//------Noting the cpu start time
	time_Ugb = clock();
	//------Initialization of the values to be determined
	for (int j = 1; j <= JMAX; j++)
	{
		for (int i = 1; i <= IMAX-1; i++)
		{
			//------Define cell size and distances----------
			double dx_p, dx_e, dx_w;
			double dy_p, dy_n, dy_s;
			dx_p	= 0.5*(dx[i] + dx[i+1]);
			dy_p	= dy[j];
			dx_e	= dx[i+1];
			dx_w	= dx[i];
			dy_n	= 0.5*(dy[j] + dy[j+1]);
			if (j == 1)	dy_s	= dy[j]; 
			else		dy_s	= 0.5*(dy[j-1] + dy[j]);	
			//------Computation of values at face-centers--
			double mue_e = mue(i+1,j);
			double mue_w = mue(i,j);
			double mue_n = 0.25*(mue(i,j) + mue(i+1,j) + mue(i,j+1) + mue(i+1,j+1));
			double mue_s = 0.25*(mue(i,j) + mue(i+1,j) + mue(i,j-1) + mue(i+1,j-1));
			//------Adjustment of values for boundaries--
			if (j == 1)				// Bottom Wall
			{
				mue_s	= 0.25*(mue(i,j) + mue(i+1,j) + mue(i,j) + mue(i+1,j));
			}
			else if (j == JMAX)		// Top Wall
			{
				mue_n 	= 0.25*(mue(i,j) + mue(i+1,j) + mue(i,j) + mue(i+1,j));
			}
			//------X Staggered velocity at cell faces-----
			double U_p	= Ug[i][j];
			double U_pe	= Ug[i+1][j];
			double U_pw	= Ug[i-1][j];
			double U_pn	= Ug[i][j+1];
			double U_ps	= Ug[i][j-1];						
			double U_e	= 0.5*(U_p + U_pe);
			double U_w	= 0.5*(U_p + U_pw);
			double U_n	= 0.5*(U_p + U_pn);
			double U_s	= 0.5*(U_p + U_ps);
			//------Y Staggered velocity at cell faces-----
			double V_nw	= Vg[i][j];
			double V_ne	= Vg[i+1][j];
			double V_sw	= Vg[i][j-1];
			double V_se	= Vg[i+1][j-1];		
			double V_n	= 0.50*(V_nw + V_ne);
			double V_s	= 0.50*(V_sw + V_se);
			double V_e	= 0.50*(V_ne + V_se);
			double V_w	= 0.50*(V_nw + V_sw);
			//------ Convective Fluxes in the X direction--	
			Fe_x[i][j]	= (1.00/dx_p) * 0.50 * rho_g * U_e;
			Fw_x[i][j]	= (1.00/dx_p) * 0.50 * rho_g * U_w;
			Fn_x[i][j]	= (1.00/dy_p) * 0.50 * rho_g * V_n;
			Fs_x[i][j]	= (1.00/dy_p) * 0.50 * rho_g * V_s;
			//------ Diffusive Fluxes in the X direction--
			De_x[i][j]	= (mue_e/(dx_p*dx_e));
			Dw_x[i][j]	= (mue_w/(dx_p*dx_w));
			Dn_x[i][j]	= (mue_n/(dy_p*dy_n));
			Ds_x[i][j]	= (mue_s/(dy_p*dy_s));
			//------ Sources and Sinks in the X direction--
			//------ TURBULENT KINETIC ENERGY COMPUTATION
			double KE_E				= KE[i+1][j];
			double KE_P 			= KE[i][j];
			double Turbulent_KE_term= (2.0/3.0) * rho_g *  (KE_E - KE_P)/dx_p; 
			double Ext_forces 		= 0.00;
			Source_x[i][j]			= Ext_forces - Turbulent_KE_term; 
			// Determination of Peclet number
			Peclet_x[i][j]	= (U_e + U_w)* dx_p/(mue_e + mue_w);
			Peclet_y[i][j]	= (V_n + V_s)* dy_p/(mue_n + mue_s);
			// Stores value of the previous iteration
			Previous[i][j]	= 0.00;						
		}
	}
	int explicit_convection_flag = 1;
	int explicit_diffusion_flag = 0;
	//------Loop for the determination of u* (f)
	while (bound>tol)
	{
		Conv_U++;
		for (int j = 1; j <= JMAX; j++)
		{
			for (int i = 1; i <= IMAX; i++)
			{
				if(i == 0)					// Left Wall
				{
					f[0][j]		=	u_left;
				}
				else if (i == IMAX)			// Right Wall
				{
					f[IMAX][j]	=	u_right;
				}
				else
				{	
					//------Define the velocity variables-----------
					double U_p	= Ug[i][j];
					double U_pn	= Ug[i][j+1];
					double U_ps	= Ug[i][j-1];
					double U_pe	= Ug[i+1][j];
					double U_pw	= Ug[i-1][j];
					//------Define the intermediate velocity variables-----------
					double f_p	= f[i][j];
					double f_pn	= f[i][j+1];
					double f_ps	= f[i][j-1];
					double f_pe	= f[i+1][j];
					double f_pw	= f[i-1][j];			
					//------Adjustment of values for boundaries--
					if (j == 1)						// Bottom Wall
					{
						f_ps	= 2.00 * u_bottom - f_p;
					}
					else if (j == JMAX)				// Top Wall
					{
						f_pn	= 2.00 * u_top - f_p;
					}
					//=============================================== FLUX COMPUTATION ==============================================//
					// Convective fluxes
					double Fe		= 	Fe_x[i][j];
					double Fw		= 	Fw_x[i][j];
					double Fn		= 	Fn_x[i][j];
					double Fs		= 	Fs_x[i][j];
					// Diffusive fluxes
					double De		= 	De_x[i][j];
					double Dw		= 	Dw_x[i][j];
					double Dn		= 	Dn_x[i][j];
					double Ds		= 	Ds_x[i][j];
					//=============================================== CONVECTION COMPUTATION ========================================//
					//------Convection term(Implicit)
					double alpha_x= 1.00;
					if (Peclet_x[i][j] < -2) 	alpha_x = -1;
					else if (Peclet_x[i][j] > 2)	alpha_x = 1;
					else							alpha_x = 0;					
					double alpha_y= 1.00;
					if (Peclet_y[i][j] < -2)		alpha_y = -1;
					else if (Peclet_y[i][j] > 2)	alpha_y = 1;
					else							alpha_y = 0;

					double DUUDX, DVUDY;
					if( explicit_convection_flag == 1)
					{
						DUUDX		= ( Fe * ((U_p + U_pe) + alpha_x*(U_p - U_pe)) - Fw * ((U_p + U_pw) + alpha_x*(U_pw - U_p)) );
						DVUDY		= ( Fn * ((U_p + U_pn) + alpha_y*(U_p - U_pn)) - Fs * ((U_p + U_ps) + alpha_y*(U_ps - U_p)) );
					}
					else
					{
						DUUDX		= ( Fe * ((f_p + f_pe) + alpha_x*(f_p - f_pe)) - Fw * ((f_p + f_pw) + alpha_x*(f_pw - f_p)) );
						DVUDY		= ( Fn * ((f_p + f_pn) + alpha_y*(f_p - f_pn)) - Fs * ((f_p + f_ps) + alpha_y*(f_ps - f_p)) );
					}							
					double Convec_term		= DUUDX + DVUDY;
					//================================================DIFFUSION COMPUTATION ========================================//
					//------Diffusion term(Implicit)
					double D2UDX2, D2UDY2;
					if( explicit_diffusion_flag == 1)
					{
						D2UDX2		= (Dw * U_pw) + (De * U_pe);
						D2UDY2		= (Ds * U_ps) + (Dn * U_pn);
					}
					else
					{
						D2UDX2		= (Dw * f_pw) + (De * f_pe);
						D2UDY2		= (Ds * f_ps) + (Dn * f_pn);
					}
					double Diffusive_term	= D2UDX2 + D2UDY2;
					//================================================SOURCE COMPUTATION =========================================//
					double Source_term		= Source_x[i][j]; 
					//================================================PREVIOUS TIME TERM COMPUTATION ============================//
					double Previous_Time_term	= (rho_g/dt1) * U_p;
					//------Computation of u*---------------------
					double aP;
					if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 1) )
					{
						aP	= 	(rho_g/dt1); 	
					}
					else if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 0) )
					{
						aP	= 	(rho_g/dt1) + (De + Dw + Dn + Ds); 
					}
					else if( ( explicit_convection_flag == 0) && ( explicit_diffusion_flag == 1) )
					{
						aP	= 	(rho_g/dt1) + (1.00 + alpha_x)*Fe - (1.00 - alpha_x)*Fw + (1.00 + alpha_y)*Fn - (1.00 - alpha_y)*Fs;
					}
					else
					{
						aP	= 	(rho_g/dt1) + (1.00 + alpha_x)*Fe - (1.00 - alpha_x)*Fw + (1.00 + alpha_y)*Fn - (1.00 - alpha_y)*Fs + (De + Dw + Dn + Ds); 
					}
					if (aP == 0)
					{
						f[i][j]	= 0.0;	// If the grid point is in the solid phase
					}
					else
					{			
						double omega_u = 1.00;
						f[i][j] = (1.0 - omega_u)*f_p + (omega_u/aP)*(Previous_Time_term - Convec_term + Diffusive_term + Source_term);
					}					
				}
			}
		}		
		//------Enforce a state-based boundary condition-----------------------------
		enforce_physics();
		//------Computation of residue to determine convergence----------------------	
		bound	=	0.00;							
		for (int i = 0;i <= IMAX; i++)
		{
			for (int j = 1;j <= JMAX; j++)
			{
				bound				= bound + fabs(	f[i][j] - Previous[i][j]);
				Previous[i][j]	= f[i][j];
			}
		}
		if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 1) ) bound = -tol;
	}
	//------Noting the cpu end time
	time_Uga = clock();
	//------Enforce a state-based boundary condition
	enforce_physics();
	//------Noting the total cpu time for this computation
	int nthreads = 1;
	cpu_time_U_gas = ((double) (time_Uga - time_Ugb)) / (nthreads*CLOCKS_PER_SEC);
}
//----------------------------------------------------------------------
void v_update()
{
	/*------Local variables declaration--------------------------------------------------------*/
	double bound = 100.0;
	double Fe_y[IMAX+1][JMAX+1], Fw_y[IMAX+1][JMAX+1];
	double Fn_y[IMAX+1][JMAX+1], Fs_y[IMAX+1][JMAX+1];
	double De_y[IMAX+1][JMAX+1], Dw_y[IMAX+1][JMAX+1];
	double Dn_y[IMAX+1][JMAX+1], Ds_y[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1];
	double Peclet_y[IMAX+1][JMAX+1];
	double Source_y[IMAX+1][JMAX+1];
	double Previous[IMAX+1][JMAX+1];
	/*------Noting the cpu start time --------------------------------------------------------*/
	time_Vgb = clock();
	/*------Initialization of the values to be determined-------------------------------------*/
	for (int j = 1; j <= JMAX-1; j++)
	{
		for (int i = 1; i <= IMAX; i++)
		{		
			//------Define cell size and distances----------
			double dx_p, dx_e, dx_w;
			double dy_p, dy_n, dy_s;
			dx_p	= dx[i];
			dy_p	= 0.5*(dy[j] + dy[j+1]);
			dx_e	= 0.5*(dx[i] + dx[i+1]);
			if (i == 1) 	dx_w	=	dx[i]; 
			else 			dx_w	=	0.5*(dx[i-1] + dx[i]);
			dy_n	= dy[j+1];
			dy_s	= dy[j] ;
			//------Computation of values at face-centers--
			double mue_e = 0.25*(mue(i,j+1) + mue(i,j) + mue(i+1,j+1) + mue(i+1,j));
			double mue_w = 0.25*(mue(i,j+1) + mue(i,j) + mue(i-1,j+1) + mue(i-1,j));
			double mue_n = mue(i,j+1);
			double mue_s = mue(i,j);
			//------Adjustment of values for boundaries--
			if (i == 1)				// Left Wall
			{
				mue_w	= 0.25*(mue(i,j+1) + mue(i,j) + mue(i,j+1) + mue(i,j));
			}
			else if(j == JMAX-1)	// Right Wall
			{
				mue_e	= 0.25*(mue(i,j+1) + mue(i,j) + mue(i,j+1) + mue(i,j));
			}		
			//------X Staggered velocity at cell faces-----
			double U_nw	= Ug[i-1][j+1];
			double U_ne	= Ug[i][j+1];
			double U_sw	= Ug[i-1][j];
			double U_se	= Ug[i][j];
			double U_e	= 0.500*(U_ne + U_se);
			double U_w	= 0.500*(U_nw + U_sw);			
			double U_n	= 0.500*(U_ne + U_nw);
			double U_s	= 0.500*(U_se + U_sw);
			//------Y Staggered velocity at cell faces-----
			double V_p	= Vg[i][j];
			double V_pe	= Vg[i+1][j];
			double V_pw	= Vg[i-1][j];
			double V_pn	= Vg[i][j+1];
			double V_ps	= Vg[i][j-1];
			double V_e	= 0.5*(V_p + V_pe);
			double V_w	= 0.5*(V_p + V_pw);			
			double V_n	= 0.5*(V_p + V_pn);
			double V_s	= 0.5*(V_p + V_ps);
			//------ Convective Fluxes in the Y direction--	
			Fe_y[i][j]	= (1.00/dx_p) * 0.50 * rho_g * U_e;
			Fw_y[i][j]	= (1.00/dx_p) * 0.50 * rho_g * U_w;
			Fn_y[i][j]	= (1.00/dy_p) * 0.50 * rho_g * V_n;
			Fs_y[i][j]	= (1.00/dy_p) * 0.50 * rho_g * V_s;
			//------ Diffusive Fluxes in the Y direction--
			De_y[i][j]	= (mue_e/(dx_p*dx_e));
			Dw_y[i][j]	= (mue_w/(dx_p*dx_w));
			Dn_y[i][j]	= (mue_n/(dy_p*dy_n));
			Ds_y[i][j]	= (mue_s/(dy_p*dy_s));
			//------ Sources and Sinks in the Y direction--
			//================================================TURBULENT KINETIC ENERGY COMPUTATION====================//
			double KE_N = KE[i][j+1];
			double KE_P = KE[i][j];
			double Turbulent_KE_term 	= (2.0/3.0) * rho_g *  (KE_N - KE_P)/dy_p;	
			double Ext_forces 		= 0.00;
			Source_y[i][j]	= Ext_forces - Turbulent_KE_term;

			// Determination of Peclet number
			Peclet_x[i][j]	= (U_e + U_w)* dx_p/(mue_e + mue_w);
			Peclet_y[i][j]	= (V_n + V_s)* dy_p/(mue_n + mue_s);
			// Stores value of the previous iteration
			Previous[i][j]	= 0.00;		
		}
	}
	int explicit_convection_flag = 1;
	int explicit_diffusion_flag = 0;
	//------Loop for the determination of v* (g)
	while (bound>tol)
	{
		Conv_V++;
		for (int j = 0; j <= JMAX; j++)
		{
			for (int i = 1; i <= IMAX; i++)
			{
				if(j == 0)					// Bottom Wall
				{
					g[i][0]		=	v_bottom;
				}
				else if (j == JMAX)			// Front Wall
				{
					g[i][JMAX]	=	v_top;
				}
				else
				{	
					//------Define the velocity variables
					double V_p	= Vg[i][j];
					double V_pn	= Vg[i][j+1];
					double V_ps	= Vg[i][j-1];
					double V_pe	= Vg[i+1][j];
					double V_pw	= Vg[i-1][j];
					//------Define the intermediate velocity variables
					double g_p	= g[i][j];
					double g_pn	= g[i][j+1];
					double g_ps	= g[i][j-1];
					double g_pe	= g[i+1][j];
					double g_pw	= g[i-1][j];
					//------Adjustment of values for boundaries
					if (i == 1)						// Left Wall
					{
						g_pw	= 2.00 * v_left - g_p;
					}
					else if(i == IMAX)				// Right Wall
					{
						g_pe	= 2.00 * v_right - g_p;
					}
					//=============================================== FLUX COMPUTATION ==============================================//
					// Convective fluxes
					double Fe		= 	Fe_y[i][j];
					double Fw		= 	Fw_y[i][j];
					double Fn		= 	Fn_y[i][j];
					double Fs		= 	Fs_y[i][j];
					// Diffusive fluxes
					double De		= 	De_y[i][j];
					double Dw		= 	Dw_y[i][j];
					double Dn		= 	Dn_y[i][j];
					double Ds		= 	Ds_y[i][j];
					//=============================================== CONVECTION COMPUTATION ========================================//
					//------Convection term(Implicit)
					double alpha_x= 1.00;
					if (Peclet_x[i][j] < -2) 	alpha_x = -1;
					else if (Peclet_x[i][j] > 2)	alpha_x = 1;
					else							alpha_x = 0;					
					double alpha_y= 1.00;
					if (Peclet_y[i][j] < -2)		alpha_y = -1;
					else if (Peclet_y[i][j] > 2)	alpha_y = 1;
					else							alpha_y = 0;

					double DUVDX, DVVDY;
					if( explicit_convection_flag == 1)
					{
						DUVDX		= ( Fe * ((V_p + V_pe) + alpha_x*(V_p - V_pe)) - Fw * ((V_p + V_pw) + alpha_x*(V_pw - V_p)) );
						DVVDY		= ( Fn * ((V_p + V_pn) + alpha_y*(V_p - V_pn)) - Fs * ((V_p + V_ps) + alpha_y*(V_ps - V_p)) );
					}
					else
					{
						DUVDX		= ( Fe * ((g_p + g_pe) + alpha_x*(g_p - g_pe)) - Fw * ((g_p + g_pw) + alpha_x*(g_pw - g_p)) );
						DVVDY		= ( Fn * ((g_p + g_pn) + alpha_y*(g_p - g_pn)) - Fs * ((g_p + g_ps) + alpha_y*(g_ps - g_p)) );
					}							
					double Convec_term		= DUVDX + DVVDY;
					//------Diffusion term(Implicit)
					double D2VDX2, D2VDY2;
					if( explicit_diffusion_flag == 1)
					{
						D2VDX2		= (Dw * V_pw) + (De * V_pe);
						D2VDY2		= (Ds * V_ps) + (Dn * V_pn);
					}
					else
					{
						D2VDX2		= (Dw * g_pw) + (De * g_pe);
						D2VDY2		= (Ds * g_ps) + (Dn * g_pn);
					}
					double Diffusive_term	= D2VDX2 + D2VDY2;
					//================================================SOURCE COMPUTATION =========================================//
					double Source_term		= Source_y[i][j]; 
					//================================================PREVIOUS TIME TERM COMPUTATION ============================//
					double Previous_Time_term	= (rho_g/dt1) * V_p;
					//------Computation of v*
					double aP;
					if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 1) )
					{
						aP	= 	(rho_g/dt1); 	
					}
					else if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 0) )
					{
						aP	= 	(rho_g/dt1) + (De + Dw + Dn + Ds); 
					}
					else if( ( explicit_convection_flag == 0) && ( explicit_diffusion_flag == 1) )
					{
						aP	= 	(rho_g/dt1) + (1.00 + alpha_x)*Fe - (1.00 - alpha_x)*Fw + (1.00 + alpha_y)*Fn - (1.00 - alpha_y)*Fs;
					}
					else
					{
						aP	= 	(rho_g/dt1) + (1.00 + alpha_x)*Fe - (1.00 - alpha_x)*Fw + (1.00 + alpha_y)*Fn - (1.00 - alpha_y)*Fs + (De + Dw + Dn + Ds); 
					}
					if (aP == 0)
					{
						g[i][j]	= 0.0;	// If the grid point is in the solid phase
					}
					else
					{			
						double omega_v = 1.00;
						g[i][j] = (1.0 - omega_v)*g_p + (omega_v/aP)*(Previous_Time_term - Convec_term + Diffusive_term + Source_term);
					}					
				}					
			}
		}
		//------Enforce a state-based boundary condition
		enforce_physics();
		//------Computation of residue to determine convergence
		bound	=	0.00;							
		for (int i = 1;i <= IMAX; i++)
		{
			for (int j = 0;j <= JMAX; j++)
			{
				bound			= bound + fabs(	g[i][j] - Previous[i][j]);
				Previous[i][j]	= g[i][j];
			}
		}
		if( ( explicit_convection_flag == 1) && ( explicit_diffusion_flag == 1) ) bound = -tol;
	}
	//------Noting the cpu end time
	time_Vga = clock();
	//------Enforce a state-based boundary condition
	enforce_physics();
	//------Noting the total cpu time for this computation
	int nthreads = 1;
	cpu_time_V_gas = ((double) (time_Vga - time_Vgb)) / (nthreads*CLOCKS_PER_SEC);
}
//----------------------------------------------------------------------
void P_update()
{
	//------Local variables declaration
	double bound = 100.0;
	double Fe_z[IMAX+1][JMAX+1], Fw_z[IMAX+1][JMAX+1];
	double Fn_z[IMAX+1][JMAX+1], Fs_z[IMAX+1][JMAX+1];

	double De_z[IMAX+1][JMAX+1], Dw_z[IMAX+1][JMAX+1];
	double Dn_z[IMAX+1][JMAX+1], Ds_z[IMAX+1][JMAX+1];

	double Peclet_x[IMAX+1][JMAX+1];
	double Peclet_y[IMAX+1][JMAX+1];
	double Previous[IMAX+1][JMAX+1];
	double time_Pb, time_Pa, cpu_time_Pressure;
	double omega_p = 1.00;
	int loop = 0;
	//------Noting the cpu start time
	time_Pb = clock();
	//------Initialization the variable that stores the previous iteration value
	for (int j = 1; j <= JMAX; j++)
	{
		for (int i = 1; i <= IMAX; i++)
		{					
			Previous[i][j] = 0.00;
		}
	}
	while (bound > tol)
	{
		Conv_P++;
		for (int j = 1; j <= JMAX; j++)
		{
			for (int i = 1; i <= IMAX; i++)
			{					
				//------Define cell size and distances
				double dx_p, dx_e, dx_w;
				double dy_p, dy_n, dy_s;
				dx_p	=	dx[i];
				dy_p	=	dy[j];
			
				if ( i == 1 ) dx_w = dx[i]; 
				else dx_w = 0.50 * (dx[i-1] + dx[i]);
				if ( i == IMAX ) dx_e = dx[i]; 
				else dx_e = 0.50 * (dx[i] + dx[i+1]);

				if ( j == 1 ) dy_s = dy[j]; 
				else dy_s = 0.50 * ( dy[j-1] + dy[j] );				
				if ( j == JMAX ) dy_n = dy[j]; 
				else dy_n = 0.50 * ( dy[j] + dy[j+1] );
				//------Define the approximate gas velocity variables
				double f_e	= f[i][j];
				double f_w	= f[i-1][j];
				double g_n	= g[i][j];
				double g_s	= g[i][j-1];
				//------Define the variables with gas pressure at previous timestep
				double P_P	= Previous[i][j];
				double P_E	= P[i+1][j];
				double P_W	= P[i-1][j];
				double P_N	= P[i][j+1]; 
				double P_S	= P[i][j-1]; 
				//------Adjustment of values for boundaries
				if(i == 1)
				{
					P_W	= P[i][j]; 
				}
				else if(i == IMAX)
				{
					P_E	= P[i][j];
				}
				
				if(j == 1)
				{
					P_S	= P[i][j]; 
				}
				else if(j == JMAX)
				{
					P_N	= P[i][j];
				}
				//------Coefficient term computation
				double aE= (1.0/(dx_p*dx_e));
				double aW= (1.0/(dx_p*dx_w));
				double aN= (1.0/(dy_p*dy_n));
				double aS= (1.0/(dy_p*dy_s));
				double aP= (aE + aW + aN + aS);
				//------Computation of pressure derivative terms
				double D2PDX2 = aW*P_W + aE*P_E;
				double D2PDY2 = aS*P_S + aN*P_N;
				//------Computation of divergence of approximate velocity(RHS)
				double DUDX	=	(f_e - f_w)/dx_p;
				double DVDY	=	(g_n - g_s)/dy_p;
				double RHS	=	(rho_g/dt1)*(DUDX + DVDY);
				//------Computation of pressure term
				P[i][j] = ((1.0 - omega_p) * P_P + (omega_p/aP) * (D2PDX2 + D2PDY2 - RHS));
			}
		}
		double max_P = -1e-6;
		bound	=	0.00;
		for (int i = 1; i <= IMAX; i++)
		{
			for (int j = 1; j <= JMAX; j++)
			{
				bound = bound + fabs( P[i][j] - Previous[i][j] );
				Previous[i][j] = P[i][j];//relaxp * p[i][j][k];
				if (P[i][j] > max_P)
				{
					max_P = P[i][j];
				}	
			}
		}	
		loop = loop + 1;
		if ((loop%10000) == 0)printf("Pressure Evolution :: %d %lf %lf %lf\n", loop, bound, max_P, omega_p);
		
		if (loop == 100000)
		{
			printf("Pressure Evolution :: Forced Exit \n");
			break;
		}
	}
	//------Convergence met. p has been determined
	//------Noting the cpu end time
	time_Pa = clock();
	//------Noting the total cpu time for this computation
	int nthreads = 1;
	cpu_time_Pressure = ((double) (time_Pa - time_Pb)) / (nthreads*CLOCKS_PER_SEC);
}
//----------------------------------------------------------------------
void correction()
{
	for (int i = 1 ;i <= IMAX; i++)
	{
		for (int j = 1 ;j <= JMAX; j++)
		{
			//------Define cell size and distances----------------------------------------------------
			double dx_p, dx_e, dx_w;
			double dy_p, dy_n, dy_s;
			dx_p	=	dx[i];
			dy_p	=	dy[j];
			
			if ( i == 1 ) dx_w = dx[i]; 
			else dx_w = 0.50 * (dx[i-1] + dx[i]);
			if ( i == IMAX ) dx_e = dx[i]; 
			else dx_e = 0.50 * (dx[i] + dx[i+1]);

			if ( j == 1 ) dy_s = dy[j]; 
			else dy_s = 0.50 * ( dy[j-1] + dy[j] );				
			if ( j == JMAX ) dy_n = dy[j]; 
			else dy_n = 0.50 * ( dy[j] + dy[j+1] );			
			//------Define approximate velocity ----------------------------------------------------------
			double f_p 	= f[i][j];
			double g_p 	= g[i][j];
			//------Define gas pressure-------------------------------------------------------------------
			double P_P 	= P[i][j];
			double P_E 	= P[i+1][j];
			double P_N 	= P[i][j+1];
			//------Compute the corrected gas X velocity--------------------------------------------------
			Ug[i][j] = f_p - (dt1/rho_g) * ( P_E - P_P )/dx_e;
			Vg[i][j] = g_p - (dt1/rho_g) * ( P_N - P_P )/dy_n;
			//------Adjustment of values for boundaries---------------------------------------------------
			Ug[0][j]		= u_left;	
			Ug[IMAX+1][j]	= u_right;					
			Vg[i][0]		= v_bottom;	
			Vg[i][JMAX+1]	= v_top;	
		}
	}
}
//----------------------------------------------------------------------
void KE_update()
{
	/*------Local variables declaration--------------------------------------------------------*/
	double bound = 100.0;
	double Fe_k[IMAX+1][JMAX+1], Fw_k[IMAX+1][JMAX+1];
	double Fn_k[IMAX+1][JMAX+1], Fs_k[IMAX+1][JMAX+1];
	double De_k[IMAX+1][JMAX+1], Dw_k[IMAX+1][JMAX+1];
	double Dn_k[IMAX+1][JMAX+1], Ds_k[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1];
	double Peclet_y[IMAX+1][JMAX+1];
	double Source_k[IMAX+1][JMAX+1];
	double Previous[IMAX+1][JMAX+1];
	double time_KEb, time_KEa, cpu_time_KE;
	double omega_ke = 1.00;
	/*------Noting the cpu start time ---------------------------------------------------------*/
	time_KEb = clock();
	/*------Initialization the variable that stores the previous iteration value---------------*/
	for (int j = 1; j <= JMAX; j++)
	{
		for (int i = 1; i <= IMAX; i++)
		{						
			//------Define cell size and distances--------------------------------
			double dx_p, dx_e, dx_w;
			double dy_p, dy_n, dy_s;
			dx_p	=	dx[i];
			dy_p	=	dy[j];
			
			if ( i == 1 ) dx_w = dx[i]; 
			else dx_w = 0.50 * (dx[i-1] + dx[i]);
			if ( i == IMAX ) dx_e = dx[i]; 
			else dx_e = 0.50 * (dx[i] + dx[i+1]);

			if ( j == 1 ) dy_s = dy[j]; 
			else dy_s = 0.50 * ( dy[j-1] + dy[j] );				
			if ( j == JMAX ) dy_n = dy[j]; 
			else dy_n = 0.50 * ( dy[j] + dy[j+1] );
			//------Computation of values at face-centers-----------------------------
			double mut_P	= mut(i,j)/sigma_k;
			double mut_N	= mut(i,j+1)/sigma_k;
			double mut_S	= mut(i,j-1)/sigma_k;
			double mut_W	= mut(i-1,j)/sigma_k;
			double mut_E	= mut(i+1,j)/sigma_k;
			
			double mut_n	= 0.50 * (mut_P + mut_N); 
			double mut_s	= 0.50 * (mut_P + mut_S);
			double mut_w	= 0.50 * (mut_P + mut_W); 
			double mut_e	= 0.50 * (mut_P + mut_E);	
			/*------Define the gas velocity variables---------------------------------*/
			double U_w	= Ug[i-1][j]; 
			double U_e	= Ug[i][j];
			double V_n	= Vg[i][j]; 
			double V_s	= Vg[i][j-1];
			//------ Convective Fluxes 
			Fe_k[i][j]	= (1.00/dx_p) * 0.50 * rho_g * U_e;
			Fw_k[i][j]	= (1.00/dx_p) * 0.50 * rho_g * U_w;
			Fn_k[i][j]	= (1.00/dy_p) * 0.50 * rho_g * V_n;
			Fs_k[i][j]	= (1.00/dy_p) * 0.50 * rho_g * V_s;
			//------ Diffusive Fluxes 
			De_k[i][j]	= (mut_e/(dx_p*dx_e));
			Dw_k[i][j]	= (mut_w/(dx_p*dx_w));
			Dn_k[i][j]	= (mut_n/(dy_p*dy_n));
			Ds_k[i][j]	= (mut_s/(dy_p*dy_s));
			//------ Sources and Sinks
			//------Turbulence generation term computation----------------------------
			double DE_P	= DE[i][j];
			
			double U_Ne	= Ug[i][j+1];
			double U_Nw	= Ug[i-1][j+1]; 
			double U_Se	= Ug[i][j-1];
			double U_Sw	= Ug[i-1][j-1]; 
			
			double V_nE	= Vg[i+1][j]; 
			double V_sE	= Vg[i+1][j-1]; 
			double V_sW	= Vg[i-1][j-1]; 
			double V_nW	= Vg[i-1][j]; 
					
			double DUDX	= (1.0/dx_p)*(U_e - U_w);
			double DVDY	= (1.0/dy_p)*(V_n - V_s);
			
			double DUDY	= (0.25/dy_p)*((U_Ne + U_Nw) - (U_Se + U_Sw));
			double DVDX	= (0.25/dx_p)*((V_nE - V_sE) - (V_sW - V_nW));
			
			double Generation_term = mut_P * ( 2.0 * (sqr(DUDX)+sqr(DVDY)) + sqr(DUDY + DVDX));
			//------Turbulent dissipation term computation----------------------------
			double Dissipation_term 	= rho_g * DE_P;	
			Source_k[i][j] = Generation_term - Dissipation_term;
			// Stores value of the previous iteration				
			Previous[i][j] = 0.00;
		}
	}
	/*------Loop for the determination of ke------------------------------------------------------------------*/
	while (bound>tol)
	{
		Conv_KE++;
		for (int j = 1; j <= JMAX; j++)
		{
			for (int i = 1; i <= IMAX; i++)
			{					
				/*------Define the variables with gas turbulent ke at previous timestep---*/
				double KE_P				= KE[i][j];	
				/*------Define the intermediate gas turbulent ke variables-----------*/
				double KE1_P			= ke1[i][j];
				double KE1_N			= ke1[i][j+1]; 
				double KE1_S			= ke1[i][j-1]; 
				double KE1_E			= ke1[i+1][j];
				double KE1_W			= ke1[i-1][j]; 
				//------Convection term using upwind (Implicit)------------
				double alpha = 1.00;
				double Fe				= Fe_k[i][j];
				double Fw				= Fw_k[i][j];
				double Fn				= Fn_k[i][j];
				double Fs				= Fs_k[i][j];
				double DUKDX			= ( Fe * ((KE1_P + KE1_E) + alpha*(KE1_P - KE1_E)) - Fw * ((KE1_W + KE1_P) + alpha*(KE1_W - KE1_P)) );
				double DVKDY			= ( Fn * ((KE1_P + KE1_N) + alpha*(KE1_P - KE1_N)) - Fs * ((KE1_S + KE1_P) + alpha*(KE1_S - KE1_P)) ); 
				double Convec_term		= DUKDX + DVKDY;
				/*------Diffusion term using upwind (Implicit for stability)--------------*/
				double De				= De_k[i][j];
				double Dw				= Dw_k[i][j];
				double Dn				= Dn_k[i][j];
				double Ds				= Ds_k[i][j];
				double D2KDX2			= (KE1_E * De) + (KE1_W * Dw);
				double D2KDY2			= (KE1_N * Dn) + (KE1_S * Ds);
				double Diffusive_term	= D2KDX2 + D2KDY2;
				//------ Sources and Sinks
				double Source_term		= Source_k[i][j];
				/*------Previous time term computation------------------------------------*/
				double Previous_Time_term		= (rho_g/dt1) * KE_P;
				/*------Computation of gas turbulent ke-----------------------------------*/
				double aP = (rho_g/dt1) + (De + Dw + Dn + Ds);
				if (aP == 0) // Added recently
				{
					ke1[i][j]	= 0.00;	// If the grid point is in the solid phase
				}
				else
				{
					ke1[i][j]	= (1.0-omega_ke)*KE_P + (omega_ke/aP)*(Previous_Time_term - Convec_term + Diffusive_term + Source_term);
				}
			}
		}
		/*------Computation of residue to determine convergence----------------------*/	
		bound	=	0.00;
		for (int i = 1; i <= IMAX ;i++)
		{
			for (int j = 1; j <= JMAX ;j++)
			{
				bound			=	bound + fabs(ke1[i][j] - Previous[i][j]);
				Previous[i][j]	=	ke1[i][j];//relaxke * ke1[i][j];
			}
		}		
	}
	/*------Convergence met. ke has been determined--------------------------------------*/		
	for (int i = 1; i <= IMAX ;i++)
	{
		for (int j = 1; j <= JMAX ;j++)
		{
			KE[i][j] = ke1[i][j];
		}
	}
	/*------Noting the cpu end time ----------------------------------------------------------------------*/
	time_KEa = clock();
	/*------Enforce a state-based boundary condition------------------------------------------------------*/
	enforce_physics();
	/*------Noting the total cpu time for this computation------------------------------------------------*/
	int nthreads = 1;
	cpu_time_KE = ((double) (time_KEa - time_KEb)) / (nthreads*CLOCKS_PER_SEC);		
}
//----------------------------------------------------------------------
void DE_update()
{
	/*------Local variables declaration--------------------------------------------------------*/
	double bound = 100.0;
	double Fe_e[IMAX+1][JMAX+1], Fw_e[IMAX+1][JMAX+1];
	double Fn_e[IMAX+1][JMAX+1], Fs_e[IMAX+1][JMAX+1];
	double De_e[IMAX+1][JMAX+1], Dw_e[IMAX+1][JMAX+1];
	double Dn_e[IMAX+1][JMAX+1], Ds_e[IMAX+1][JMAX+1];
	double Peclet_x[IMAX+1][JMAX+1];
	double Peclet_y[IMAX+1][JMAX+1];
	double Peclet_z[IMAX+1][JMAX+1];
	double Source_e[IMAX+1][JMAX+1];
	double Previous[IMAX+1][JMAX+1];
	double time_DEb, time_DEa, cpu_time_DE;
	double omega_de = 1.00;
	/*------Noting the cpu start time ---------------------------------------------------------*/
	time_DEb = clock();
	/*------Initialization the variable that stores the previous iteration value---------------*/
	for (int j = 1; j <= JMAX; j++)
	{
		for (int i = 1; i <= IMAX; i++)
		{						
			//------Define cell size and distances--------------------------------
			double dx_p, dx_e, dx_w;
			double dy_p, dy_n, dy_s;
			dx_p	=	dx[i];
			dy_p	=	dy[j];
			
			if ( i == 1 ) dx_w = dx[i]; 
			else dx_w = 0.50 * (dx[i-1] + dx[i]);
			if ( i == IMAX ) dx_e = dx[i]; 
			else dx_e = 0.50 * (dx[i] + dx[i+1]);

			if ( j == 1 ) dy_s = dy[j]; 
			else dy_s = 0.50 * ( dy[j-1] + dy[j] );				
			if ( j == JMAX ) dy_n = dy[j]; 
			else dy_n = 0.50 * ( dy[j] + dy[j+1] );
			//------Computation of values at face-centers-----------------------------
			double mut_P	= mut(i,j)/sigma_e;
			double mut_N	= mut(i,j+1)/sigma_e;
			double mut_S	= mut(i,j-1)/sigma_e;
			double mut_W	= mut(i-1,j)/sigma_e;
			double mut_E	= mut(i+1,j)/sigma_e;
		
			double mut_n	= 0.50 * (mut_P + mut_N); 
			double mut_s	= 0.50 * (mut_P + mut_S);
			double mut_w	= 0.50 * (mut_P + mut_W); 
			double mut_e	= 0.50 * (mut_P + mut_E);	
			//------Define the gas velocity variables---------------------------------
			double U_w	= Ug[i-1][j];
			double U_e	= Ug[i][j];
			double V_n	= Vg[i][j]; 
			double V_s	= Vg[i][j-1];
			//------ Convective Fluxes 
			Fe_e[i][j]	= (1.00/dx_p) * 0.50 * rho_g * U_e;
			Fw_e[i][j]	= (1.00/dx_p) * 0.50 * rho_g * U_w;
			Fn_e[i][j]	= (1.00/dy_p) * 0.50 * rho_g * V_n;
			Fs_e[i][j]	= (1.00/dy_p) * 0.50 * rho_g * V_s;
			//------ Diffusive Fluxes 
			De_e[i][j]	= (mut_e/(dx_p*dx_e));
			Dw_e[i][j]	= (mut_w/(dx_p*dx_w));
			Dn_e[i][j]	= (mut_n/(dy_p*dy_n));
			Ds_e[i][j]	= (mut_s/(dy_p*dy_s));
			//------ Sources and Sinks
			//------Turbulence generation term computation----------------------------
			double KE_P	= KE[i][j];
			double DE_P	= DE[i][j];
			
			double U_Ne	= Ug[i][j+1];
			double U_Nw	= Ug[i-1][j+1]; 
			double U_Se	= Ug[i][j-1];
			double U_Sw	= Ug[i-1][j-1]; 
			
			double V_nE	= Vg[i+1][j]; 
			double V_sE	= Vg[i+1][j-1]; 
			double V_sW	= Vg[i-1][j-1]; 
			double V_nW	= Vg[i-1][j]; 
				
			double DUDX	= (1.0/dx_p)*(U_e - U_w);
			double DVDY	= (1.0/dy_p)*(V_n - V_s);
			
			double DUDY	= (0.25/dy_p)*((U_Ne + U_Nw) - (U_Se + U_Sw));
			double DVDX	= (0.25/dx_p)*((V_nE - V_sE) - (V_sW - V_nW));
					
			double Generation_term = (DE_P/KE_P) * c_1 * (mut_P * ( 2.0 * (sqr(DUDX)+sqr(DVDY)) + sqr(DUDY + DVDX)));
			//------Turbulent dissipation term computation----------------------------
			double Dissipation_term 	= (DE_P/KE_P) * c_2 * rho_g * DE_P;	
			Source_e[i][j] = Generation_term - Dissipation_term;
			// Stores value of the previous iteration				
			Previous[i][j] = 0.00;
		}
	}
	//------Loop for the determination of ke------------------------------------------------------------------
	while (bound>tol)
	{
		Conv_DE++;
		for (int j = 1; j <= JMAX; j++)
		{
			for (int i = 1; i <= IMAX; i++)
			{					
				/*------Define the variables with gas turbulent ke at previous timestep---*/
				double DE_P				= DE[i][j];	
				/*------Define the intermediate gas turbulent ke variables-----------*/
				double DE1_P			= de1[i][j];
				double DE1_N			= de1[i][j+1]; 
				double DE1_S			= de1[i][j-1]; 
				double DE1_E			= de1[i+1][j];
				double DE1_W			= de1[i-1][j]; 
				//------Convection term using upwind (Implicit)------------
				double alpha = 1.00;
				double Fe				= Fe_e[i][j];
				double Fw				= Fw_e[i][j];
				double Fn				= Fn_e[i][j];
				double Fs				= Fs_e[i][j];
				double DUeDX			= ( Fe * ((DE1_P + DE1_E) + alpha*(DE1_P - DE1_E)) - Fw * ((DE1_W + DE1_P) + alpha*(DE1_W - DE1_P)) );
				double DVeDY			= ( Fn * ((DE1_P + DE1_N) + alpha*(DE1_P - DE1_N)) - Fs * ((DE1_S + DE1_P) + alpha*(DE1_S - DE1_P)) ); 
				double Convec_term		= DUeDX + DVeDY;
				/*------Diffusion term using upwind (Implicit for stability)--------------*/
				double De				= De_e[i][j];
				double Dw				= Dw_e[i][j];
				double Dn				= Dn_e[i][j];
				double Ds				= Ds_e[i][j];
				double D2eDX2			= (DE1_E * De) + (DE1_W * Dw);
				double D2eDY2			= (DE1_N * Dn) + (DE1_S * Ds);
				double Diffusive_term	= D2eDX2 + D2eDY2;
				//------ Sources and Sinks
				double Source_term		= Source_e[i][j];
				/*------Previous time term computation------------------------------------*/
				double Previous_Time_term		= (rho_g/dt1) * DE_P;
				/*------Computation of gas turbulent ke-----------------------------------*/
				double aP = (rho_g/dt1) + (De + Dw + Dn + Ds);
				if (aP == 0) // Added recently
				{
					de1[i][j]	= 0.00;	// If the grid point is in the solid phase
				}
				else
				{
					de1[i][j]	= (1.0-omega_de)*DE_P + (omega_de/aP)*(Previous_Time_term - Convec_term + Diffusive_term + Source_term);
				}
			}
		}
		/*------Computation of residue to determine convergence----------------------*/	
		bound	=	0.00;
		for (int i = 1; i <= IMAX ;i++)
		{
			for (int j = 1; j <= JMAX ;j++)
			{
				bound				=	bound + fabs(de1[i][j] - Previous[i][j]);
				Previous[i][j]	=	de1[i][j];//relaxde * de1[i][j][k];
			}
		}		
	}
	/*------Convergence met. ke has been determined--------------------------------------*/		
	
	/*------??--------------*/
	for (int i = 1; i <= IMAX ;i++)
	{
		for (int j = 1; j <= JMAX ;j++)
		{
			DE[i][j] = de1[i][j];
		}
	}
	/*------Noting the cpu end time ----------------------------------------------------------------------*/
	time_DEa = clock();
	/*------Enforce a state-based boundary condition------------------------------------------------------*/
	enforce_physics();
	/*------Noting the total cpu time for this computation------------------------------------------------*/
	int nthreads = 1;
	cpu_time_DE = ((double) (time_DEa - time_DEb)) / (nthreads*CLOCKS_PER_SEC);		
}
//====================================================================================================================
void Steady_State_report()
{
	FILE *hp;
	//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS
	hp=fopen("2.3.Hydrodynamics.dat","w");
	fprintf(hp,"TITLE = Hydrodynamic state of the domain\n");
	fprintf(hp, "VARIABLES=\"X (m)\",\"Y (m)\",\"U gas (m/s)\",\"V gas (m/s)\",\"P (Pa)\",\"KE (J)\",\"DE (J)\",\"State\"\n");
	fprintf(hp, "ZONE T = \"%lf\", I=%d, J=%d F=POINT\n",solution_time, IMAX, JMAX);	
	
		for(int j = 1; j <= JMAX; j++)
		{
			for(int i = 1; i <= IMAX;i++)
			{
				fprintf(hp,"%e %e %e %e %e %e %e %d\n", x[i], y[j], Ug[i][j], Vg[i][j], P[i][j], KE[i][j], DE[i][j], State[i][j]);	
			}
		}
	fclose(hp);
}
/*
//====================================================================================================================
// User-defined Header files
#include "A.2.3D_U_update.h"
#include "A.3.3D_V_update.h"
#include "A.4.3D_W_update.h"
#include "A.5.3D_P_update.h"
#include "A.6.3D_KE_update.h"
#include "A.7.3D_DE_update.h"
//====================================================================================================================
void grid_generation(void);
void domain_conditions(void);
void initialize(void);
void enforce_physics(void);
void correction(void);
void reports(void);
void Steady_State_report(void);
void Transient_report(int);
//====================================================================================================================
int main ()
{
	// Domain Information
	Lx = 5.0;
	Ly = 50.0;
	Lz = 5.0;
	IMAX = 10;
	JMAX = 50;
	KMAX = 10;
	// Gas Properties 
	rho_g	= 1.178;
	mu_g	= 1.983e-5;
	// Turbulence Model
	c_mu	= 0.09;
	c_1		= 1.44;
	c_2		= 1.92;
	sigma_k	= 1.00;
	sigma_e = 1.22;
	// Properties at the left wall
	u_left	= 0.00;
	v_left	= 0.00;
	w_left	= 0.00;
	// Properties at the right wall
	u_right	= 0.00;
	v_right	= 0.00;
	w_right	= 0.00;
	// Properties at the bottom wall
	u_bottom= 0.00;
	v_bottom= 10.00;
	w_bottom= 0.00;
	// Properties at the top wall
	u_top	= 0.00;
	v_top	= 0.00;
	w_top	= 0.00;
	// Properties at the back wall
	u_back	= 0.00;
	v_back	= 0.00;
	w_back	= 0.00;
	// Properties at the front wall
	u_front	= 0.00;
	v_front	= 0.00;
	w_front	= 0.00;
	// Time step
	dt1 = 1e-6;
	
	// Grid Generation
	grid_generation();
	// Initial and Boundary Conditions
	domain_conditions();
	// Initialize Hydrodynamics
	initialize();
	enforce_physics();
	
	for (int nstep = 0; nstep < 10000; nstep++)
	{
		if( (nstep%100) == 0 ){
			printf("Solution time  :: %lf\n", dt1*nstep);
			Transient_report(nstep);
		}
		u_update();
		v_update();
		w_update();
		P_update();
		correction();
		enforce_physics();
		
		KE_update();
		DE_update();
		enforce_physics();
		Steady_State_report();
	}
	// Reports
	reports();
	return 0;
}
//====================================================================================================================
void grid_generation()
{
	int i, j, k;	
	printf("Grid Resolution :: IMAX=%d x JMAX=%d x KMAX=%d \n",IMAX, JMAX, KMAX);
	printf("Grid Type :: UNIFORM CARTESIAN GRID\n");
	printf("---------------------------------------------------\n\n");
	//------------x-direction---------------------------
	x[0] = - 0.50 * Lx/(IMAX);
	for (int i = 1; i <= IMAX; i++)
	{
		dx[i] = Lx/(IMAX);
		if(i == 1)
			x[i] = 0.50 * dx[i];
		else
			x[i] = x[i-1] + 0.50 * (dx[i] + dx[i-1]);
	} 
	x[IMAX+1] = Lx + 0.50 * Lx/(IMAX);	
	//------------y-direction---------------------------
	y[0]= - 0.50 * Ly/(JMAX);
	for (int j = 1; j <= JMAX; j++)
	{
		dy[j] = Ly/(JMAX);
		if(j == 1) 
			y[j] = 0.50 * dy[j];
		else 
			y[j] = y[j-1] + 0.50 * (dy[j] + dy[j-1]);
	}
	y[JMAX+1] = Ly + 0.50 * Ly/(JMAX);  
	//------------z-direction---------------------------
	z[0]= - 0.50 * Lz/(KMAX);
	for (int k = 1; k <= KMAX; k++)
	{
		dz[k] = Lz/(KMAX); 
		if(k == 1) 
			z[k] = 0.50 * dz[k];
		else 
			z[k] = z[k-1] + 0.50 * (dz[k] + dz[k-1]);
	}
	z[KMAX+1] = Lz + 0.50 * Lz/(KMAX);
}
//====================================================================================================================
void domain_conditions()
{
	//------------Defining Initial State----------------
	for(int i = 0; i<=IMAX+1;i++)
	{
		for(int j = 0; j<=JMAX+1; j++)
		{
			for(int k = 0; k<=KMAX+1; k++)
			{	
				State[i][j][k] = 0;
			}
		}
	}
	//------------Defining Boundary Walls----------------
	for(int j = 0; j<=JMAX+1; j++)
	{
		for(int k = 0; k<=KMAX+1; k++)
		{		
			if (State[0][j][k] == 0)		State[0][j][k]		= -1;					//Left Surface
			if (State[IMAX+1][j][k] == 0)	State[IMAX+1][j][k] = -1;					//Right Surface
		}
	}
	
	for(int i = 0; i<=IMAX+1;i++)
	{			
		for(int k = 0; k<=KMAX+1; k++)
		{		
				if(State[i][0][k] == 0)		State[i][0][k]		= -1;					//Bottom Surface
				if(State[i][JMAX+1][k] == 0)State[i][JMAX+1][k] = -1;					//Top Surface
		}
	}

	for(int i = 0; i<=IMAX+1;i++)
	{				
		for(int j = 0; j<=JMAX+1; j++)
		{
			if(State[i][j][0] == 0)			State[i][j][0]		= 1;					//Back Surface
			if(State[i][j][KMAX+1] == 0)	State[i][j][KMAX+1] = 2;					//Front Surface
		}
	}
	
	/*
	double CC_x = 1.00;
	double CC_y = 2.50;	
	double CC_z = 5.00;	
	
	double XL_brick = 2.00;
	double YL_brick = 1.00;
	double ZL_brick = 1.00;
		
	int ii1	= (CC_x - 0.5*XL_brick)/dx[1];
	int ii2	= (CC_x + 0.5*XL_brick)/dx[1];
	
	int jj1	= (CC_y - 0.5*YL_brick)/dy[1];
	int jj2	= (CC_y + 0.5*YL_brick)/dy[1];
	
	int kk1	= (CC_z - 0.5*ZL_brick)/dz[1];
	int kk2	= (CC_z + 0.5*ZL_brick)/dz[1];
	
	for(int i = 0; i < IMAX+1; i++)
	{
		for(int j = 0; j < JMAX+1; j++)
		{
			for(int k = 0;k < KMAX+1; k++)
			{
				if ( (i >= ii1 && i <= ii2)&& (j >= jj1 && j <= jj2) && (k >= kk1 && k <= kk2) )
				{
					State[i][j][k]= 1;
				}
			}
		}
	}
	
}
//====================================================================================================================
void initialize()
{
	double vin = 10.0;
	double lm = 0.005;
	for(int i = 0; i<=IMAX+1; i++)
	{
		for(int j = 0; j<=JMAX+1; j++)
		{
			for(int k = 0; k<=KMAX+1; k++)
			{	
		
				Ug[i][j][k]		= 0.00;
				Vg[i][j][k]		= 0.00;
				Wg[i][j][k]		= 0.00;
				P[i][j][k]		= 0.00;
				KE[i][j][k]		= 0.009*sqr(vin);
				DE[i][j][k]		= pow(0.009*sqr(vin),1.5)/lm;
				epfs[i][j][k]	= 0.00;
				epfd[i][j][k]	= 0.00;

				vfrac[i][j][k]	= 0.40;
				
				f[i][j][k]		= Ug[i][j][k];
				g[i][j][k]		= Vg[i][j][k];
				h[i][j][k]		= Wg[i][j][k];				
				ke1[i][j][k]	= KE[i][j][k];
				de1[i][j][k]	= DE[i][j][k];
			}
		}
	}
}
//====================================================================================================================
void enforce_physics()
{
	for(int i = 0; i <= IMAX+1; i++)
	{
		for(int j = 0; j <= JMAX+1; j++)
		{
			for(int k = 0; k <= KMAX+1; k++)
			{	
				if (State[i][j][k] == -1)			// WALL REGIONS
				{				
					Ug[i][j][k] = 0.00;
					Vg[i][j][k] = 0.00;
					Wg[i][j][k] = 0.00;	
					
					f[i][j][k] = Ug[i][j][k];
					g[i][j][k] = Vg[i][j][k];
					h[i][j][k] = Wg[i][j][k];
				}
				if (State[i][j][k] == 1)			// INLET REGIONS
				{				
					Ug[i][j][k] = 0.00;
					Vg[i][j][k] = 0.00;
					Wg[i][j][k] = 0.00;
					
					f[i][j][k] = Ug[i][j][k];
					g[i][j][k] = Vg[i][j][k];
					h[i][j][k] = Wg[i][j][k];
				}
				else if (State[i][j][k] == 2)		// OUTLET REGIONS
				{				
					P[i][j][k] = 0.00;
				}
				
				// X direction boundaries
				if(i == 0)					// Left Wall
				{
					Ug[0][j][k]		=	u_left;
					f[0][j][k]		=	u_left;
				}
				else if (i == IMAX+1)			// Right Wall
				{
					Ug[IMAX+1][j][k]=	u_right;
					f[IMAX+1][j][k]	=	u_right;
				}
				// Y direction boundaries
				if(j == 0)					// Bottom Wall
				{
					Vg[i][0][k]		=	v_bottom;
					g[i][0][k]		=	v_bottom;
				}
				else if (j == JMAX+1)			// Top Wall
				{
					Vg[i][JMAX+1][k]=	v_top;
					g[i][JMAX+1][k]	=	v_top;
				}
				// Z direction boundaries
				if(k == 0)					// Back Wall
				{
					Wg[i][j][0]		=	w_back;
					h[i][j][0]		=	w_back;
				}
				else if (k == KMAX+1)			// Front Wall
				{
					Wg[i][j][KMAX+1]=	w_front;
					h[i][j][KMAX+1]	=	w_front;
				}
			}
		}
	}
}
//====================================================================================================================
void correction()
{
	for (int i = 1;i <= IMAX; i++)
	{
		for (int j = 1;j <= JMAX; j++)
		{
			for (int k = 1;k <= KMAX; k++)
			{
				//------Define cell size and distances----------------------------------------------------
				double dx_p, dx_e, dx_w;
				double dy_p, dy_n, dy_s;
				double dz_p, dz_f, dz_b;
				dx_p	=	dx[i];
				dy_p	=	dy[j];
				dz_p	=	dz[k];
				
				if ( i == 1 ) dx_w = dx[i]; 
				else dx_w = 0.50 * (dx[i-1] + dx[i]);
				if ( i == IMAX ) dx_e = dx[i]; 
				else dx_e = 0.50 * (dx[i] + dx[i+1]);

				if ( j == 1 ) dy_s = dy[j]; 
				else dy_s = 0.50 * ( dy[j-1] + dy[j] );				
				if ( j == JMAX ) dy_n = dy[j]; 
				else dy_n = 0.50 * ( dy[j] + dy[j+1] );
				
				if ( k == 1 ) dz_b = dz[k]; 
				else dz_b = 0.50 * ( dz[k-1] + dz[k] );
				if ( k == KMAX ) dz_f = dz[k]; 
				else dz_f = 0.50 * ( dz[k] + dz[k+1] );
				//------Define approximate velocity ----------------------------------------------------------
				double f_p 	= f[i][j][k];
				double g_p 	= g[i][j][k];
				double h_p 	= h[i][j][k];
				//------Define gas pressure-------------------------------------------------------------------
				double P_P 	= P[i][j][k];
				double P_E 	= P[i+1][j][k];
				double P_N 	= P[i][j+1][k];
				double P_F 	= P[i][j][k+1];
				//------Compute the corrected gas X velocity--------------------------------------------------
				Ug[i][j][k] = f_p - (dt1/rho_g) * ( P_E - P_P )/dx_e;
				Vg[i][j][k] = g_p - (dt1/rho_g) * ( P_N - P_P )/dy_n;
				Wg[i][j][k] = h_p - (dt1/rho_g) * ( P_F - P_P )/dz_f;
				//------Adjustment of values for boundaries---------------------------------------------------
				Ug[0][j][k]			= u_left;	
				Ug[IMAX+1][j][k]	= u_right;					
				Vg[i][0][k]			= v_bottom;	
				Vg[i][JMAX+1][k]	= v_top;	
				Wg[i][j][0]			= w_back;	
				Wg[i][j][KMAX+1]	= w_front;	
			}
		}
	}
}
//====================================================================================================================
void reports()
{
	FILE *fp;
	fp = fopen("2.1.Grid Definition.dat","w");
	
		fprintf(fp, "VARIABLES=\"X\",\"Y\",\"Z\"\n");
		fprintf(fp, "ZONE T = \"Cell Centered for Scalars\"\n");	
		for(int i = 0; i<=IMAX+1;i++)
		{
			for(int j = 0; j<=JMAX+1; j++)
			{
				for(int k = 0; k<=KMAX+1; k++)
				{
					double x_p	= x[i];
					double y_p	= y[j];
					double z_p	= z[k];
					fprintf(fp,"%e %e %e\n",x_p, y_p, z_p);
				}
			}
		}
		
		fprintf(fp, "ZONE T = \"Staggered in X direction\"\n");	
		for(int i = 0; i<=IMAX; i++)
		{
			for(int j = 0; j<=JMAX+1; j++)
			{
				for(int k = 0; k<=KMAX+1; k++)
				{
					double x_p;					
					if(i == 0)
					{
						double dx_p	= 0.5*(dx[i] + dx[i+1]);
						x_p	= x[i] +  dx_p;
					}
					else
					{
						x_p	= x[i] +  0.5 * dx[i];
					}
					double y_p	= y[j];
					double z_p	= z[k];
					fprintf(fp,"%e %e %e\n",x_p,y_p,z_p);
				}
			}
		}
		
		fprintf(fp, "ZONE T = \"Staggered in Y direction\"\n");	
		for(int i = 0; i<=IMAX+1;i++)
		{
			for(int j = 0; j<=JMAX; j++)
			{
				for(int k = 0; k<=KMAX+1; k++)
				{
					double y_p;
					
					double x_p	= x[i];
					if(j == 0)
					{
						double dy_p	= 0.5*(dy[j] + dy[j+1]);	
						y_p	= y[j] + dy_p;
					}
					else
					{
						y_p	= y[j] + 0.50 * dy[j];
					}
					
					double z_p	= z[k];
					fprintf(fp,"%e %e %e\n",x_p, y_p, z_p);
				}
			}
		}
		
		fprintf(fp, "ZONE T = \"Staggered in Z direction\"\n");	
		for(int i = 0; i<=IMAX+1;i++)
		{
			for(int j = 0; j<=JMAX+1; j++)
			{
				for(int k = 0; k<=KMAX; k++)
				{
					double z_p;
					double x_p	= x[i];
					double y_p	= y[j];
					if(k == 0)
					{
						double dz_p	= 0.5*(dz[k] + dz[k+1]);
						z_p	= z[k] + dz_p;
					}
					else
					{
						z_p	= z[k] + 0.50*dz[k];
					}					
					fprintf(fp,"%e %e %e\n",x_p, y_p, z_p);
				}
			}
		}	
	fclose(fp);
	
	FILE *gp;
	//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS
	double solution_time = 0.00;
	gp=fopen("2.2.Domain_State.dat","w");
	fprintf(gp,"TITLE = State of the domain\n");
	fprintf(gp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Z (m)\",\"State\"\n");
	fprintf(gp, "ZONE T = \"%lf\", I=%d, J=%d, K=%d F=POINT\n",solution_time, IMAX+2, JMAX+2, KMAX+2);	
	
		for(int k = 0; k <= KMAX+1; k++)
		{
			for(int j = 0; j <= JMAX+1; j++)
			{
				for(int i = 0; i <= IMAX+1;i++)
				{
					fprintf(gp,"%e %e %e %d\n", x[i], y[j], z[k], State[i][j][k]);	
				}
			}
		}
	fclose(gp);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Steady_State_report()
{
	FILE *hp;
	//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS
	hp=fopen("2.3.Hydrodynamics.dat","w");
	fprintf(hp,"TITLE = Hydrodynamic state of the domain\n");
	fprintf(hp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Z (m)\",\"U gas (m/s)\",\"V gas (m/s)\",\"W gas (m/s)\",\"P (Pa)\",\"KE (J)\",\"DE (J)\",\"State\"\n");
	fprintf(hp, "ZONE T = \"%lf\", I=%d, J=%d, K=%d F=POINT\n",solution_time, IMAX, JMAX, KMAX);	
	
		for(int k = 1; k <= KMAX; k++)
		{
			for(int j = 1; j <= JMAX; j++)
			{
				for(int i = 1; i <= IMAX;i++)
				{
					fprintf(hp,"%e %e %e %e %e %e %e %e %e %d\n", x[i], y[j], z[k], Ug[i][j][k], Vg[i][j][k], Wg[i][j][k], P[i][j][k], KE[i][j][k], DE[i][j][k], State[i][j][k]);	
				}
			}
		}
	fclose(hp);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Transient_report(int N)
{
	double solution_time = N*dt1;
	FILE *hp;
	//-------------------------HYDRODYNAMIC PROPERTIES AT INTERMITTENT INTERVALS
	hp=fopen("2.4.Transient_Hydrodynamics.dat","a");
	if(N == 0)
	{
		fprintf(hp,"TITLE = Hydrodynamic state of the domain\n");
		fprintf(hp, "VARIABLES=\"X (m)\",\"Y (m)\",\"Z (m)\",\"U gas (m/s)\",\"V gas (m/s)\",\"W gas (m/s)\",\"P (Pa)\",\"KE (J)\",\"DE (J)\",\"State\"\n");
	}
	fprintf(hp, "ZONE T = \"%lf\", I=%d, J=%d, K=%d F=POINT\n",solution_time, IMAX, JMAX, KMAX);	
		for(int k = 1; k <= KMAX; k++)
		{
			for(int j = 1; j <= JMAX; j++)
			{
				for(int i = 1; i <= IMAX; i++)
				{
					fprintf(hp,"%e %e %e %e %e %e %e %e %e %d\n", x[i], y[j], z[k], Ug[i][j][k], Vg[i][j][k], Wg[i][j][k], P[i][j][k], KE[i][j][k], DE[i][j][k], State[i][j][k]);	
				}
			}
		}
	fclose(hp);
}

*/
//========================FUNCTIONS=============================================================================================
double sqr(double x)
{ 
  return (x*x); 
} 
// Viscosity
double mut(int i, int j)
{
   return ( c_mu*rho_g*(sqr(KE[i][j])/DE[i][j]) );
}
double mue(int i, int j)
{
   return ( mu_g + mut(i,j) );
}