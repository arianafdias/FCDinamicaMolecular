#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <numeric>
#include <iomanip> 
#include <string> 
#include <stdlib.h> 
using namespace std;

//global vars
double G = 4.49235*pow(10,-15); // pc Msolar^-1 (pc/ano)^2
double deltaT = 2*pow(10,4); // "erro numerico"   // 4.1 -> 2*pow(10,4)
double pcConv = 1.02269032*pow(10, -6); //conversão em parsec
struct vetor {double x; double y; double z;};
struct estrela {vetor x0; vetor x1; vetor x2; vetor v; double m;};



double distancia(vetor estrela1, vetor estrela2){

	double d = sqrt(fabs(pow((estrela1.x - estrela2.x),2.0) + pow((estrela1.y - estrela2.y),2.0) + pow((estrela1.z - estrela2.z),2.0)));
	
	return d;
}


double potencial( int npart, vector<estrela> sistema, int t){

 	double Ep = 0;
	for(int i = 0; i < npart; i++ ){
		for(int j = i+1; j < npart; j++){
			if (t == 0){
			Ep += (G*sistema[i].m*sistema[j].m)/(distancia(sistema[i].x0, sistema[j].x0));
			} else {
			Ep += (G*sistema[i].m*sistema[j].m)/(distancia(sistema[i].x1, sistema[j].x1));
			}
		}
	}
	
	return -Ep;
}

double cinetica(int npart, vector<estrela> sistema){

	// ---- VELOCIDADE DO CENTRO DE MASSA ---- //
	double mTotal = 0;
	for(int i = 0; i < npart; i++){
		mTotal += sistema[i].m;
		}
	double sumx = 0;
	double sumy = 0;
	double sumz = 0;
	for(int i = 0; i < npart; i++){
		sumx += sistema[i].m*sistema[i].v.x;
		sumy += sistema[i].m*sistema[i].v.y;
		sumz += sistema[i].m*sistema[i].v.z;
	}
	vector<double> vCM = {sumx/mTotal,sumy/mTotal,sumz/mTotal};
		
	
	// ---- ENERGIA CINÉTICA ---- //
	double somatorio = 0;
	for (int i = 0; i < npart; i++){
		double vi = pow(sistema[i].v.x,2) + pow(sistema[i].v.y,2) + pow(sistema[i].v.z,2);
		somatorio += 0.5 * sistema[i].m * vi;
	}
	
	double Ek = 0.5 * mTotal * (pow(vCM[0],2)+pow(vCM[1],2)+pow(vCM[2],2)) + somatorio;
	
	return Ek;
	
}


vetor calculo_aceleracao(int npart, vector<estrela> sistema, int t, int particulaEmQuestao){

	vetor a;
	a.x = 0;
	a.y = 0;
	a.z = 0;
	vetor F;
	F.x = 0;
	F.y = 0;
	F.z = 0;
	
	for (int i = 0; i < npart; i++){
		if (i == particulaEmQuestao) continue;
		
		if (t == 1) {
			F.x += (G * (sistema[particulaEmQuestao].m * sistema[i].m)* (sistema[particulaEmQuestao].x0.x - sistema[i].x0.x)) / pow(distancia(sistema[i].x0, sistema[particulaEmQuestao].x0), 3.0);
			F.y += (G * (sistema[particulaEmQuestao].m * sistema[i].m)* (sistema[particulaEmQuestao].x0.y - sistema[i].x0.y)) / pow(distancia(sistema[i].x0, sistema[particulaEmQuestao].x0), 3.0);
			F.z += (G * (sistema[particulaEmQuestao].m * sistema[i].m)* (sistema[particulaEmQuestao].x0.z - sistema[i].x0.z)) / pow(distancia(sistema[i].x0, sistema[particulaEmQuestao].x0), 3.0);
		} 
		else{
			F.x += (G * (sistema[particulaEmQuestao].m * sistema[i].m)* (sistema[particulaEmQuestao].x1.x - sistema[i].x1.x)) / pow(distancia(sistema[i].x1, sistema[particulaEmQuestao].x1), 3.0);
			F.y += (G * (sistema[particulaEmQuestao].m * sistema[i].m)* (sistema[particulaEmQuestao].x1.y - sistema[i].x1.y)) / pow(distancia(sistema[i].x1, sistema[particulaEmQuestao].x1), 3.0);
			F.z += (G * (sistema[particulaEmQuestao].m * sistema[i].m)* (sistema[particulaEmQuestao].x1.z - sistema[i].x1.z)) / pow(distancia(sistema[i].x1, sistema[particulaEmQuestao].x1), 3.0);
		}
		
		
		
		
	}
	a.x = -F.x / sistema[particulaEmQuestao].m;
	a.y = -F.y / sistema[particulaEmQuestao].m;
	a.z = -F.z / sistema[particulaEmQuestao].m;	
	
	
	
	return a;
}


vetor calculo_x1 ( int npart, vector<estrela> sistema, int t, int particulaEmQuestao){
	
	vetor x1;
	vetor a = calculo_aceleracao(npart, sistema, t, particulaEmQuestao);
	x1.x = sistema[particulaEmQuestao].x0.x + sistema[particulaEmQuestao].v.x * deltaT + 0.5 * a.x * pow(deltaT,2);
	x1.y = sistema[particulaEmQuestao].x0.y + sistema[particulaEmQuestao].v.y * deltaT + 0.5 * a.y * pow(deltaT,2);
	x1.z = sistema[particulaEmQuestao].x0.z + sistema[particulaEmQuestao].v.z * deltaT + 0.5 * a.z * pow(deltaT,2);
	
	return x1;
}


vetor calculo_xn ( int npart, vector<estrela> sistema, int t, int particulaEmQuestao){
	
	vetor xn;
	vetor a = calculo_aceleracao(npart, sistema, t, particulaEmQuestao);
	
	xn.x = 2.0*sistema[particulaEmQuestao].x1.x - sistema[particulaEmQuestao].x0.x + a.x*pow(deltaT,2.0);
	xn.y = 2.0*sistema[particulaEmQuestao].x1.y - sistema[particulaEmQuestao].x0.y + a.y*pow(deltaT,2.0);
	xn.z = 2.0*sistema[particulaEmQuestao].x1.z - sistema[particulaEmQuestao].x0.z + a.z*pow(deltaT, 2.0); 
	
	return xn;
}


vetor velocity(vector<estrela> sistema, int particulaEmQuestao){
	
	vetor V;
	V.x = (sistema[particulaEmQuestao].x2.x - sistema[particulaEmQuestao].x0.x) /(2.0*deltaT); 
	V.y = (sistema[particulaEmQuestao].x2.y - sistema[particulaEmQuestao].x0.y) /(2.0*deltaT);
	V.z = (sistema[particulaEmQuestao].x2.z - sistema[particulaEmQuestao].x0.z) /(2.0*deltaT); 
	
	return V;	
}


// ---------------- EXERCICIO 4.1 ---------------- //

int ex41(){

	int npart = 2; //massas das estrelas em Msolar
	vector<double> EnergiaC;
	vector<double> EnergiaP;
	
	// ------ DEFINIÇÃO DO SISTEMA ------ //
	vector<estrela> sistema;
	
	estrela new_star1;
	new_star1.x0.x = -1;
	new_star1.x0.y = 0;
	new_star1.x0.z = 0;
	
	new_star1.x1.x = 0;
	new_star1.x1.y = 0;
	new_star1.x1.z = 0;
	
	new_star1.x2.x = 0;
	new_star1.x2.y = 0;
	new_star1.x2.z = 0;
	
	new_star1.v.x = 0;
	new_star1.v.y = 2*pow(10,-2)*pcConv;
	new_star1.v.z = 0;
	
	new_star1.m = 1;
	sistema.push_back(new_star1);
	
	
	
	estrela new_star2;
	new_star2.x0.x = 1;
	new_star2.x0.y = 0;
	new_star2.x0.z = 0;
	
	new_star2.x1.x = 0;
	new_star2.x1.y = 0;
	new_star2.x1.z = 0;
	
	new_star2.x2.x = 0;
	new_star2.x2.y = 0;
	new_star2.x2.z = 0;
	
	new_star2.v.x = 0;
	new_star2.v.y = -2*pow(10,-2)*pcConv;
	new_star2.v.z = 0;
	
	new_star2.m = 1;
	sistema.push_back(new_star2);
	
	
	// SISTEMA SITUA-SE EM T = 0
	// PODE-SE CALCULAR AS ENERGIAS CINÉTICA E POTENCIAL
	int currentT = 0;
	EnergiaP.push_back(potencial(npart, sistema, currentT));
	EnergiaC.push_back(cinetica(npart, sistema));
	
	
	
	currentT ++;
	// SISTEMA SITUA-SE EM T = 1
	// CALCULA-SE AS POSIÇOES EM T = 1
	for(int n = 0; n < npart; n++){
		sistema[n].x1 = calculo_x1(npart, sistema, currentT, n);
	}
	
	
	
	currentT ++;
	// SISTEMA SITUA-SE EM T = 2
	// CALCULA-SE AS POSIÇOES EM T = 2
	for(int n = 0; n < npart; n++){
		sistema[n].x2 = calculo_xn(npart, sistema, currentT, n);
		sistema[n].v = velocity(sistema, n);
	}
	// PODE-SE CALCULAR AS ENERGIAS CINÉTICA E POTENCIAL
	EnergiaP.push_back(potencial(npart, sistema, currentT));
	EnergiaC.push_back(cinetica(npart, sistema));
	
	ofstream movfile("movV2.txt");
	for(int k = 0; k < pow(10,4); k++ ){
		currentT++;
		for( int n = 0; n < npart; n++){
			sistema[n].x0 = sistema[n].x1;
			sistema[n].x1 = sistema[n].x2;
		}	
		
		for(int n = 0; n < npart; n++){
			sistema[n].x2 = calculo_xn(npart, sistema, currentT, n);
			sistema[n].v = velocity(sistema, n);
		}
		movfile << std::fixed << std::setprecision(4) << sistema[0].x2.x << "\t" << sistema[0].x2.y << "\t" << sistema[1].x2.x << "\t" << sistema[1].x2.y << endl; 
		
	// PODE-SE CALCULAR AS ENERGIAS CINÉTICA E POTENCIAL
	EnergiaP.push_back(potencial(npart, sistema, currentT));
	EnergiaC.push_back(cinetica(npart, sistema));
	}
	
	movfile.close();
	ofstream file("Energia.txt");
	for(int i = 0; i < EnergiaC.size(); i++){
		file << EnergiaC[i] << "\t" << EnergiaP[i] << "\t" <<  EnergiaC[i] + EnergiaP[i] << endl;
	}
		
	file.close();
	
	return 0;
}

// ---------------- EXERCICIO 4.2 ---------------- //
vetor centroMassaAglomerado (vector<estrela> sistema, int t){

	vetor CM;
	
	double mTotal = 0;
	for(int i = 0; i < sistema.size(); i++){
		mTotal += sistema[i].m;
		}
		
	double sumx = 0;
	double sumy = 0;
	double sumz = 0;
	for(int i = 0; i < sistema.size(); i++){
		if (t == 0){
			sumx += sistema[i].m*sistema[i].x0.x;
			sumy += sistema[i].m*sistema[i].x0.y;
			sumz += sistema[i].m*sistema[i].x0.z;
		}else if(t == 1){
			sumx += sistema[i].m*sistema[i].x1.x;
			sumy += sistema[i].m*sistema[i].x1.y;
			sumz += sistema[i].m*sistema[i].x1.z;
		}else{
			sumx += sistema[i].m*sistema[i].x2.x;
			sumy += sistema[i].m*sistema[i].x2.y;
			sumz += sistema[i].m*sistema[i].x2.z;
		}
	}
	
	CM.x = sumx/mTotal;
	CM.y = sumy/mTotal;
	CM.z = sumz/mTotal;
			
	return CM;

}

int nrEstrelasDMenorRt (vector<estrela> sistema, int t){
	
	int nrEstrelas = 0;
	double d;
	double rT = 5; //pc
	
	vetor rCM = centroMassaAglomerado(sistema, t);
	for(int i = 0; i < sistema.size(); i++){
		if(t == 0) d = distancia(rCM,sistema[i].x0);
		else if(t == 1) d = distancia(rCM,sistema[i].x1);
		else d = distancia(rCM,sistema[i].x1);
		
		if (d < rT) nrEstrelas++;
	}

	return nrEstrelas;
}


int ex42(){

	// ---- TIRAR INFO DO FILE ---- //
	ifstream myfile;
	myfile.open("500.dat");
	int npart = 500;
	estrela newStar;
	vector<estrela> sistema(npart, newStar);
	int idx = 0;
	while (!myfile.eof()){
		myfile >> sistema[idx].x0.x;
		myfile >> sistema[idx].x0.y;
		myfile >> sistema[idx].x0.z;
		myfile >> sistema[idx].v.x;
		myfile >> sistema[idx].v.y;
		myfile >> sistema[idx].v.z;
		myfile >> sistema[idx].m;
		idx++;
	}
	
	
	
	// ------- SIMULAÇÃO ------- //
	int tMAX = 10*pow(10,6)/deltaT;
	vector<int> nrTotalEstrelasMenorRt;
	
	
	// SISTEMA SITUA-SE EM T = 0
	int currentT = 0;
	// AS POSIÇÕES INICIAIS SÃO RETIRADAS DA LEITURA DO FICHEIRO
	
	nrTotalEstrelasMenorRt.push_back(nrEstrelasDMenorRt(sistema, currentT));
	currentT++;
	
	
	// SISTEMA SITUA-SE EM T = 1
	// CALCULA-SE AS POSIÇOES EM T = 1
	for(int n = 0; n < npart; n++){
		sistema[n].x1 = calculo_x1(npart, sistema, currentT, n);
	}
	
	// CALCULA-SE O CENTRO DE MASSA DO SISTEMA EM T = 1
	// VERIFICA-SE AS DISTANCIAS DAS ESTRELAS AO CM DO AGLOMERADO
	nrTotalEstrelasMenorRt.push_back(nrEstrelasDMenorRt(sistema, currentT));
	currentT++;
	
	
	// SISTEMA SITUA-SE EM T = 2
	// CALCULA-SE AS POSIÇOES EM T = 2
	for(int n = 0; n < npart; n++){
		sistema[n].x2 = calculo_xn(npart, sistema, currentT, n);
	}
	
	// CALCULA-SE O CENTRO DE MASSA DO SISTEMA EM T = 2
	// VERIFICA-SE AS DISTANCIAS DAS ESTRELAS AO CM DO AGLOMERADO
	nrTotalEstrelasMenorRt.push_back(nrEstrelasDMenorRt(sistema, currentT));
	currentT++;
	
	
	
	// SISTEMA SITUA-SE EM T > 2
	for(int t = currentT; t < tMAX; t ++){
		for( int n = 0; n < npart; n++){
			sistema[n].x0 = sistema[n].x1;
			sistema[n].x1 = sistema[n].x2;
		}	
		
		for(int n = 0; n < npart; n++){
			sistema[n].x2 = calculo_xn(npart, sistema, t, n);
			/*
			ofstream myFile("EstrelaNr" + to_string(n) + ".txt", fstream::app);
			myFile << sistema[n].x2.x << "\t" << sistema[n].x2.y << "\t" << sistema[n].x2.z << endl;
			myFile.close();
			*/
			
			
		}
		
		nrTotalEstrelasMenorRt.push_back(nrEstrelasDMenorRt(sistema, t));
	 }
	
	ofstream file("Distruibuicoes"+to_string(npart)+".txt");
	for(int i = 0; i < nrTotalEstrelasMenorRt.size(); i++){
		file << (double) nrTotalEstrelasMenorRt[i] / (double) npart << endl;
	}
	
	
	
	file.close();



	return 0;
}


// ---------------- EXERCICIO 4.3 ---------------- //

double gerarMassa() {
	
	double m;
	double c = 0.205827;
	double random = drand48();
	double alpha1 = 0.3;
	double alpha2 = 1.3;
	double alpha3 = 2.3;
	
	if (0 < random && random <= 0.0384783) {
		m = pow((((1-alpha1)/c)*random + pow(0.01, (1-alpha1))),(1/(1-alpha1)));
	} else if (0.0384783 < random && random <= 0.610547) {
		m = pow(((1-alpha2)*( random/c - pow(0.08,(1-alpha1))/(1-alpha1) + pow(0.01,(1-alpha1))/(1-alpha1) + pow(0.08,(1-alpha2))/(1-alpha2))),(1/(1-alpha2)));
	} else if (0.610547 < random && random <= 1){
		m = pow(    (1-alpha3)*( random/ c -pow(0.08,(1-alpha1))/(1-alpha1) +pow(0.01,(1-alpha1))/(1-alpha1) -pow(0.5,(1-alpha2))/(1-alpha2) + pow(0.08,(1-alpha2))/(1-alpha2) + pow(0.5,(1-alpha3))/(1-alpha3)    ), (1 / (1-alpha3))) ;
	}
	
	return m;
}

vector<estrela> gerarPosicoesIniciais(vector<estrela> sistema){
	
	vector<double> amostrasAceitaveisParaR;
	
	while (amostrasAceitaveisParaR.size() < npart){
		double randomX = drand48();	
		double randomY = drand48();
		double OrdenadaYdeRandomX = PI(randomX);
		
		if (randomY <= OrdenadaYdeRandomX) amostrasAceitaveisParaR.push_back(randomX);
		
	}
	
	vector<double> thetas;
	vector<double> phis;
	for (int i = 0; i < npart; i++){
		double theta = 2* M_PI * amostrasAceitaveisParaR;
		thetas.push_back(theta);
		double phi = acos(2 * amostrasAceitaveisParaR - 1);
		phis.push_back(phi);
		
	}	
	
	for(int i = 0; i < npart; i++){
		sistema[i].x0.x = amostrasAceitaveisParaR[i] * cos(thetas[i]) * sin(phis[i]);
		sistema[i].x0.y = amostrasAceitaveisParaR[i] * sin(thetas[i]) * sin(phis[i]);
		sistema[i].x0.z = amostrasAceitaveisParaR[i] * cos(phis[i]);
			
	return sistema;
}



double PI(double r){
	
	double PIr = (asinh(r) - r / sqrt(26))/ ((sqrt(26) * asinh(5) - 5 )/sqrt(26)) ;
	return PIr;
	}

int ex43(){ 

	int npart = 100;
	estrela newStar;
	vector<estrela> sistema(npart, newStar);
	
	for(int i = 0; i < npart; i++){
		
		sistema[i].m = gerarMassa();
		//cout << sistema[i].m  << endl;
	
	}
	
	sistema = gerarPosicoesIniciais(sistema);
	
	
	
	
	
	
	
	
	
	
	
	
	return 0;
}

int main(){
	//ex41();
	//ex42();
	ex43();
	return 0;
}
