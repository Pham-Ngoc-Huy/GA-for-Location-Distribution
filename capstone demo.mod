/*********************************************
 * OPL 12.10.0.0 Model
 * Author: nghai
 * Creation Date: Nov 25, 2022 at 4:37:04 PM
 *********************************************/
//parameters
int H = ...; //Set of suppliers
range rangeh = 1..H; 
int I = ...; //Set of distribution centers
range rangei = 1..I; 
int J = ...; //Set of customers
range rangej = 1..J;
int K = ...; //Set of fresh food types
range rangek = 1..K;
int U[rangeh] = ...; //Maximum supply capacity of supplier  to the enterprise
int N[rangei] = ...; //Maximum inventory capacity of distribution center 
float R[rangek] = ...; //Corruption rate coefficient of fresh food , constant
float l1[rangeh][rangei] = ...; //Distance between supplier  and distribution center 
float l2[rangei][rangej] = ...; //Distance between distribution center  and customer 
int v = ...; //Average speed of vehicle
float c1[rangeh][rangei] =...; //Transportation cost per unit from supplier  to distribution center
float c2[rangei][rangej] = ...; //Transportation cost per unit from distribution center  to customer 
float A[rangei] =...; //Fixed cost of distribution center 
float B[rangei] = ...; //Operating cost of distribution center 
float f[rangei][rangej] =...; //Penalty costs of vehicle arriving late to customer 
int w=...; //The largest number of new distribution centers available to rent
float r[rangek] =...; //The average selling price of the variety  of fresh goods
float t[rangei][rangej] =...; //Time of vehicle from distribution center  to customer 
float T[rangej] = ...; //The serve time that customer  requires from the delivery center
float Tmax[rangej] = ...; //The shortest waiting time that customer  is not dissatisfied
float Tmin[rangej] = ...; //The longest waiting time that customer  is dissatisfied
float a[rangej][rangek] = ...; //The freshness satisfaction threshold for customer  of variety  of fresh goods
float b[rangej] =...;  //The time satisfaction threshold for customer .
int M = 1000000;
int gamma = 2;
float freshness[rangej][rangek];
float A1 = 0.5;
float A2 = 0.5;
float q;
float omega;
float a2 = ...; // mean (1-R[rangek])
float zmax = 24062.290125; // maximize
range x = 0.. 24062.290125; 
//decision variables
dvar int+ d1[rangeh][rangei];
dvar int+ d2[rangei][rangej];
dvar int+ d[rangej];
dvar boolean y[rangei];
dvar boolean z1[rangeh][rangei];
dvar boolean z2[rangei][rangej];
dvar boolean z3[rangej][rangek];
dvar boolean r1[rangei][rangej];
dvar boolean s[rangei][rangej];
dvar float+ s9;
dvar float+ s1; //fixed costs and operating expenses
dvar float+ s2; //transportation costs
dvar float+ s3; //delayed penalty costs
dvar float+ s10;
dvar float+ s4; //fresh food corruption costs
dvar float+ s11;
dvar float+ mean[rangej]; // the distribution center service time 
dvar float+ Y;
dvar float+ t2[rangej][rangek]; // the freshness of fresh goods
//for the modified linearize function


// Objective functions:
minimize s1 + s2 + s3 + s4;
maximize ;
maximize s11;


subject to{
//s1
s1 == sum(i in rangei)(A[i]*y[i]) + sum(i in rangei)(B[i]*y[i]);
//s2
s2 == sum(h in rangeh, i in rangei)(c1[h][i]*l1[h][i]*d1[h][i]) + sum( i in rangei, j in rangej)(c2[i][j]*l2[i][j]*d2[i][j]);
//s3
forall(i in rangei, j in rangej){
T[j] <= t[i][j] + M*(1- s[i][j]);
t[i][j] <= M*s[i][j];
}
sum(i in rangei, j in rangej) s[i][j] == 1;
forall(i in rangei, j in rangej){
 s3 == s[i][j] * f[i][j]*(t[i][j] - T[j]);
}
//s4
forall(k in rangek){
  (1-R[k])/sum(k in rangek)R[k] == a2;
  }
 
  

//objective function 10

  
//ti[i][k]
s9 == sum(h in rangeh, i in rangei, j in rangej)((A1*((l1[h][i]/v) + t[i][j]) + A2*(l1[h][i] + l2[i][j]))*z2[i][j]);
q == sum(h in rangeh, i in rangei, j in rangej)((A1*((l1[h][i]/v) + t[i][j]) + A2*(l1[h][i] + l2[i][j])));


//mean[j]
forall(i in rangei, j in rangej){
Tmin[j] <= t[i][j] + M*(1 - r1[i][j]);
t[i][j] <= Tmax[j] + M*(1 - r1[i][j]);
t[i][j] <= M*r1[i][j];
}
sum(i in rangei, j in rangej)(r1[i][j]) == 1;
sum(i in rangei, j in rangej)(r1[i][j] * ((Tmax[j]- t[i][j])/(Tmax[j]-Tmin[j]))^gamma) == Y; 

forall(j in rangej){
 mean[j] == sum(i in rangei, j in rangej)(z2[i][j] *Y);
 s11 == mean[j];
} 

// constraint 1 
forall(h in rangeh){
  sum(i in rangei) d1[h][i] <= U[h];
 }
// constraint 2 & 3 
forall(i in rangei){
  sum(h in rangeh) d1[h][i] <= N[i]*y[i];
  sum(j in rangej) d2[i][j] <= N[i]*y[i];
 }
 // constraint 4&5
 forall(j in rangej){
sum( i in rangei) y[i] <= w;
sum( i in rangei) z2[i][j] ==1;
}
// constraint 6
forall(j in rangej){
  sum( i in rangei) d2[i][j] >= d[j];
}
// constraint 7&8
forall( i in rangei, j in rangej, h in rangeh){
  z2[i][j] <= y[i];
  z1[h][i] <= y[i];
}
// constraint 9&10
forall( i in rangei, j in rangej, h in rangeh){
  d1[h][i] >= 0;
  d2[i][j] >= 0;  
}
// constraint 11
forall( i in rangei, j in rangej, k in rangek ){
  t[i][j] >= a[j][k];
}
// constraint 12


// constraint 13
forall(j in rangej){
 mean[j] >= b[j];
}
} 
//meta-heuristics -FOA: