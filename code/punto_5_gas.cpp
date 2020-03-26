#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

using namespace std;

const double Lx=60,Ly=60;
const int Nx=5,Ny=5, N=Nx*Ny;
const double pared[4]={0,Lx,0,2*Ly};

const double Zi=0.1786178958448091e0;
const double Lambda=0.2123418310626054*(-1);
const double Xi=0.06626458266981849*(-1);

const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=(1-2*(Xi+Zi));

bool choque_x[N],choque_y[N];

class Cuerpo{
	private:

		vector3D r,V,F; double m,R;

	public:

		void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0){
			r.cargue(x0,y0,0); V.cargue(Vx0,Vy0,0); m=m0; R=R0;
		}
		double fuerza_pared(double position_on_axis, double wall_position, bool *choque){
			double h=R-abs(wall_position-position_on_axis);

			if(h>=0){*choque=true; return 1e4*pow(h,1.5);}
			else if(h<0){*choque=false; return 0;}
		}
		void AgregueFuerza(vector3D dF, int counter){
			vector3D aux_Fx, aux_Fy; double F_x=0.0, F_y=0.0;

			for(int i=0; i<2; i++){
				F_x+=pow(-1.0,i)*fuerza_pared(r.x(), pared[i], &choque_x[counter]);
				F_y+=pow(-1.0,i)*fuerza_pared(r.y(), pared[i+2], &choque_y[counter]);
			}
			aux_Fx.cargue(F_x,0,0);aux_Fy.cargue(0,F_y,0);
			//cout<<F_x<<'\t'<<F_y<<endl;
			F += dF+aux_Fx+aux_Fy;
		}
		void Mueva_r(double dt,double Coeficiente){r+=V*(dt*Coeficiente);}
		void Mueva_V(double dt,double Coeficiente){V+=F*(dt*Coeficiente/m);}
		void Dibujese(void){cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";}
		void BorreFuerza(void){F.cargue(0,0,0);}
		double Getx(void){return r.x();} //Inline
		double Gety(void){return r.y();} //Inline
		double Getz(void){return r.z();} //Inline
		double GetVx(void){return V.x();} //Inline
		double GetVy(void){return V.y();} //Inline
		double GetFx(void){return F.x();} //Inline
		double GetFy(void){return F.y();} //Inline
		friend class Colisionador;
};

class Colisionador{
	public:
		void CalculeFuerzaEntre(Cuerpo & Molecula1,Cuerpo & Molecula2, int counter_i, 
			int counter_j){
			vector3D dr=Molecula2.r-Molecula1.r; double Norma2_dr=norma2(dr);

			// aux=(12*epsilon*pow(r_0, 6)*pow(Norma2_dr,-4))    *    ((pow(r_0,6)*pow(Norma2_dr, -3))-1);
			double aux=(12.0e6*pow(Norma2_dr,-4))    *    (((1.0e6)*pow(Norma2_dr, -3))-1);
			
			vector3D F2=dr*aux;
			Molecula2.AgregueFuerza(F2, counter_i);	Molecula1.AgregueFuerza(F2*(-1), counter_j);
		};
		void CalculeTodasLasFuerzas(Cuerpo * Molecula){
			int i,j;
			for(i=0;i<(N+4);i++) Molecula[i].BorreFuerza();
			for(i=0;i<N;i++)
				for(j=i+1;j<(N+4);j++)
					CalculeFuerzaEntre(Molecula[i],Molecula[j], i, j);
		};
		void mover_todo_segun_PEFRL(Cuerpo * Molecula, double dt){
			int i;
			for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Zi);
			CalculeTodasLasFuerzas(Molecula);
			for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Coeficiente1);
			for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Xi);
			CalculeTodasLasFuerzas(Molecula);
			for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Lambda);
			for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Coeficiente2);
			CalculeTodasLasFuerzas(Molecula);
			for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Lambda);
			for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Xi);
			CalculeTodasLasFuerzas(Molecula);
			for(i=0;i<N;i++) Molecula[i].Mueva_V(dt,Coeficiente1);
			for(i=0;i<N;i++) Molecula[i].Mueva_r(dt,Zi);
		};
};
//------------------ Funciones Globales -----------------
void InicieAnimacion(void);
void InicieCuadro(void);
void TermineCuadro(void);

int main(void){
	Cuerpo Molecula[N];
	Colisionador Newton;
	Crandom ran64(1);

	double t,tdibujo,dt=1.0e-3;
	int i,j;
	
	InicieAnimacion();

	//-------------------------INICIALIZACIÓN-------------------------
	double m0=1.0,    R0=2.5,    k_T=100;
	double dx=Lx/(Nx+1),    dy=Ly/(Ny+1),    x0, y0, theta, Vx0, Vy0;
	double V0=sqrt(2.0*k_T/m0);
		
	for(i=0;i<Nx;i++)
		for(j=0;j<Ny;j++){
			x0=(i+1)*dx; y0=(j+1)*dy; 
			theta=2.0*M_PI*ran64.r(); Vx0=V0*cos(theta); Vy0=V0*sin(theta);
			//---------------------(x0,y0,Vx0,Vy0, m0, R0)
			Molecula[i*Ny+j].Inicie(x0, y0, Vx0, Vy0, m0, R0);
		}
	
	for(i=0;i<N;i++){choque_x[i]=false; choque_y[i]=false;};

	//-------------------------CORRER LA SIMULACIÓN-------------------------
	double T_sim=30.0,T_eq=100.0; //t_eq=100
	/*
	double y_prom, suma, UnoSobreN=1.0/N;
	double Vx_init[N], Vy_init[N], Intensidad, delta_V;
	int num_datos=(T_sim+dt)*N/dt, Cont= 0;
	double Vel_x[num_datos];
	cout<<num_datos<<endl;
	*/
	for(t=tdibujo=0.0;t<T_eq+T_sim;t+=dt,tdibujo+=dt){
		
		if(tdibujo>0.05){
			InicieCuadro();
			for(i=0;i<N;i++) Molecula[i].Dibujese();
			TermineCuadro();
			tdibujo=0.0;
		}
		
		/*
		//--- Búsqueda de y promedio ---
		y_prom=0.0;suma=0.0;
		 
		for(i=0;i<N;i++){suma+=Molecula[i].Gety();};
		
		y_prom=suma*UnoSobreN;
		cout<<t<<'\t'<<y_prom<<endl;
		
		
		//---Velocidades en x----- 
		if(t>T_eq){
			for(i=0;i<N;i++){
				//Vel_x[Cont+i]=Molecula[i].GetVx();
				cout<<Molecula[i].GetVx()<<endl;				
				//Cont+=1;
			}
		}
		
		//--- Intensidad Choques ---
		
		for(i=0;i<N;i++){Vx_init[i]=0.0; Vy_init[i]=0.0;}

		if(t>T_eq){
			Newton.CalculeTodasLasFuerzas(Molecula);

			for(i=0;i<N;i++){
				//if(choque_x[i]) Vx_init[i]=abs(Molecula[i].GetVx());
				//if(choque_y[i]) Vy_init[i]=abs(Molecula[i].GetVy());
				if(choque_x[i]) delta_V+=abs(Molecula[i].GetFx()*dt);
				if(choque_y[i]) delta_V+=abs(Molecula[i].GetFy()*dt);
				cout<<choque_x[i]<<'\t'<<choque_y[i]<<endl;	
			}
		}
		*/
		//for(i=0;i<N;i++){choque_x[i]=false; choque_y[i]=false;};
		
		//Moverlo Segun PEFRL Orden 4
		Newton.mover_todo_segun_PEFRL(Molecula, dt);
		/*
		if(t>T_eq){
			for(i=0;i<N;i++){
				//if(choque_x[i]) delta_V+=abs(Molecula[i].GetVx()-Vx_init[i]);
				//if(choque_y[i]) delta_V+=abs(Molecula[i].GetVy()-Vy_init[i]);
			}
		cout <<'\n'<< endl;
		cout <<'\n'<< endl;
		}*/

	}
	/*
	Intensidad+=delta_V/T_sim;
	//cout << Intensidad/60.0 << endl;
	// kT=2, P= 
	// kT=4, P= 
	// kT=6, P=
	// kT=8, P=
	// kT=10, P=
	// kT=15, P=
	// kT=20, P=
	*/
	return 0;
}

//------------------ Funciones Globales -----------------
void InicieAnimacion(void){
	//cout<<"set terminal gif animate"<<endl; 
	//cout<<"set output 'pelicula.gif'"<<endl;
	cout<<"unset key"<<endl;
	cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
	cout<<"set yrange[-10:"<<2*Ly+10<<"]"<<endl;
	cout<<"set size ratio -1"<<endl;
	cout<<"set parametric"<<endl;
	cout<<"set trange [0:7]"<<endl;
	cout<<"set isosamples 12"<<endl;	
}
void InicieCuadro(void){
		cout<<"plot 0,0 ";
		cout<<" , "<<Lx/7<<"*t,0";			  //pared de abajo
		cout<<" , "<<Lx/7<<"*t,"<<2*Ly;	  //pared de arriba
		cout<<" , 0,"<<2*Ly/7<<"*t";				//pared de la izquierda
		cout<<" , "<<Lx<<","<<2*Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
		cout<<endl;
}
