#include"itensor/all.h"
#include <fstream>

using namespace std;
using namespace itensor;

int main(int argc, char *argv[])
{
  using S1 = MixedSiteSet<ElectronSite,SpinHalfSite>;
  // want a chain of lenght L, choosen small to compare with ED values
  auto L=160;
  auto N=2*L;
  auto tc=1.0;
  auto JC=0.5;
  auto U=20.0;
  auto JS=0.5;
  auto JK=0.9; 
  auto JCS=0.0;    
     
    
   printf("%d     %.3f    %.3f    %.3f    %.3f     %.3f\n",L,JC,U,JS,JK,JCS);

    //generating the mixed site set with conserved fermion number, but not bosons

    //uto sites = S1(N,{"ConserveQNs=",true});
    //auto sites = S1(N,{"ConserveNf",false, "ConserveSz",true});

     S1 sites;
     readFromFile("L=160sites",sites);
     MPS psi(sites);
     readFromFile("L=160psi",psi);

  auto ampo = AutoMPO(sites);

 // generating the Hamiltonian, making sure local fermionic operators only act on site
 // 1, 3, 5 etc
 // and itinerant fermionic  on 2, 4 ..

   
    for(int j=1; j<= N-3; j+=2)
        { 
          ampo += -tc,"Cdagup",j,"Cup",j+2;
          ampo += -tc,"Cdagup",j+2,"Cup",j;
          ampo += -tc,"Cdagdn",j,"Cdn",j+2;
          ampo += -tc,"Cdagdn",j+2,"Cdn",j;   
        }
   for(int j=1; j<= N-3; j+=2)
        {
          ampo += JC/2,"S+",j,"S-",j+2;
          ampo += JC/2,"S-",j,"S+",j+2;
          ampo +=   JC,"Sz",j,"Sz",j+2;
        }

    for(int j=1; j<= N-3; j+=2)
        {
          ampo += JS/2,"LS+",j+1,"LS-",j+3;
          ampo += JS/2,"LS-",j+1,"LS+",j+3;
          ampo +=   JS,"LSz",j+1,"LSz",j+3;
        }   
    
    for(int j = 1; j <= N-1; j += 2)
      {
        ampo += JK/2,"S+",j,"LS-",j+1;
        ampo += JK/2,"S-",j,"LS+",j+1;
        ampo +=   JK,"Sz",j,"LSz",j+1;
      }

    for(int j = 1; j <= N-3; j += 2)
      {
        ampo += JCS/2,"S+",j,"LS-",j+3;
        ampo += JCS/2,"S-",j,"LS+",j+3;
        ampo +=   JCS,"Sz",j,"LSz",j+3;

        ampo += JCS/2,"LS+",j+1,"S-",j+2;
        ampo += JCS/2,"LS-",j+1,"S+",j+2;
        ampo +=   JCS,"LSz",j+1,"Sz",j+2;
      }     

    for(int j=1; j<= N-1; j+=2)
        {
          ampo += U,"Nup",j,"Ndn",j;        
        }

    auto H=toMPO(ampo);       

    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
     /*auto state = InitState(sites);
      for(int i = 1; i <= N; i++)//local
        {
        
        if((i == 151)||(i == 153)||(i == 155)||(i == 157)||(i == 159)||(i == 161)||(i == 163)||(i == 165)||(i == 167)||(i == 169)||(i == 171)||(i == 173)) 
           {  
             println("itinerant Emp ", i); //12 hole
             state.set(i,"Emp");
           }
        else if(i%4 == 1) 
             {
              println("itinerant Up ", i);
              state.set(i,"Up");
             }
        else if(i%4 == 3)         
            { 
             println("itinerant Dn ", i);              
             state.set(i,"Dn");
            }
       
        else if(i%4 == 2)
           {
             println("spin LDn ", i); 
             state.set(i,"LDn");
           }
        else if(i%4 == 0)          
          {
             println("spin LUp ", i);     
             state.set(i,"LUp");
           }
             
        }
    auto psi = MPS(state);
    print(totalQN(psi));
    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.5f", inner(psi,H,psi) );

    auto sweeps = Sweeps(500);
    sweeps.maxdim() = 4000;
    sweeps.cutoff() = 1E-9;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-2,1E-3,1E-4,1E-5,1E-6,1E-7;
    println(sweeps);


    auto energy = dmrg(psi,H,sweeps,"Quiet");
    printfln("\nGround State Energy = %.10f",energy);
    printfln("\nUsing inner = %.10f", inner(psi,H,psi) );

    writeToFile("L=160sites",sites);
    writeToFile("L=160psi",psi);*/


/*
//LSz
ofstream SaveFile1("loc_LSz.dat");
printfln("loc_LSz= ");
  for (int i=2;i<N;i=i+2)
     {
       psi.position(i);
       auto psidag = dag(psi);                 
       auto C = elt(psi(i)*op(sites,"LSz",i)*prime(psidag(i),"Site"));
       SaveFile1<<i<<"     "<<C<<endl; 
       
      }
  SaveFile1.close();
  println("\n");

//Sz
ofstream SaveFile2("iti_Sz.dat");
printfln("iti_Sz= ");
  for (int i=1;i<=N;i=i+2)
     {
       psi.position(i);
       auto psidag = dag(psi);                 
       auto C = elt(psi(i)*op(sites,"Sz",i)*prime(psidag(i),"Site"));
       SaveFile2<<i<<"     "<<C<<endl; 
       
      }
  SaveFile2.close();
  println("\n");

//Nhole
ofstream SaveFile3("iti_hole.dat");
printfln("iti_hole= ");
  for (int i=1;i<=N;i=i+2)
     {
       psi.position(i);
       auto psidag = dag(psi);                 
       auto C = elt(psi(i)*op(sites,"Ntothole",i)*prime(psidag(i),"Site"));
       SaveFile3<<i<<"     "<<C<<endl; 
       
      }
  SaveFile3.close();
  println("\n");

//iti_ss
ofstream SaveFile5("iti_ss.dat");
printfln("iti_ss= ");    

 for (int i = 1; i <= N; i+=2)
   {
        
   for (int j = i+2; j <= N; j+=2)
      {
        
        
        psi.position(i);
        auto psidag = dag(psi);                 
        auto ir = commonIndex(psi(i), psi(i + 1), "Link");

        auto C1 = psi(i);
        C1 *= op(sites, "Sz", i);
        C1 *= prime(prime(psidag(i), "Site"), ir);            
      
        auto C2 = psi(i);
        C2 *= op(sites, "S+", i);
        C2 *= prime(prime(psidag(i), "Site"), ir);            
      

        auto C3 = psi(i);
        C3 *= op(sites, "S-", i);
        C3 *= prime(prime(psidag(i), "Site"), ir);            

            for (int k = i + 1; k < j; ++k)
            {
                 C1 *= psi(k);
                 C1 *= prime(psidag(k), "Link"); 
                 C2 *= psi(k);
                 C2 *= prime(psidag(k), "Link"); 
                 C3 *= psi(k);
                 C3 *= prime(psidag(k), "Link"); 
             }

            auto jl = commonIndex(psi(j), psi(j - 1), "Link");

            C1 *= psi(j);
            C1 *= op(sites, "Sz", j);
            C1 *= prime(prime(psidag(j), "Site"), jl);
             
            C2 *= psi(j);
            C2 *= op(sites, "S-", j);
            C2 *= prime(prime(psidag(j), "Site"), jl);
           
            C3 *= psi(j);
            C3 *= op(sites, "S+", j);
            C3 *= prime(prime(psidag(j), "Site"), jl);
            
            auto sij = elt(C1)+ elt(C2)+ elt(C3);
             
             SaveFile5<<i<<"     "<<j<<"     "<<sij<<endl;
             SaveFile5<<j<<"     "<<i<<"     "<<sij<<endl;
             
        }
    }
  SaveFile5.close();
  println("\n");

//loc_SS
ofstream SaveFile10("loc_SS.dat");
printfln("loc_SS= ");    

 for (int i = 2; i <= N; i+=2)
   {
        
   for (int j = i+2; j <= N; j+=2)
      {
        
        
        psi.position(i);
        auto psidag = dag(psi);                 
        auto ir = commonIndex(psi(i), psi(i + 1), "Link");

        auto C1 = psi(i);
        C1 *= op(sites, "LSz", i);
        C1 *= prime(prime(psidag(i), "Site"), ir);            
      
        auto C2 = psi(i);
        C2 *= op(sites, "LS+", i);
        C2 *= prime(prime(psidag(i), "Site"), ir);            
      

        auto C3 = psi(i);
        C3 *= op(sites, "LS-", i);
        C3 *= prime(prime(psidag(i), "Site"), ir);            

            for (int k = i + 1; k < j; ++k)
            {
                 C1 *= psi(k);
                 C1 *= prime(psidag(k), "Link"); 
                 C2 *= psi(k);
                 C2 *= prime(psidag(k), "Link"); 
                 C3 *= psi(k);
                 C3 *= prime(psidag(k), "Link"); 
             }

            auto jl = commonIndex(psi(j), psi(j - 1), "Link");

            C1 *= psi(j);
            C1 *= op(sites, "LSz", j);
            C1 *= prime(prime(psidag(j), "Site"), jl);
             
            C2 *= psi(j);
            C2 *= op(sites, "LS-", j);
            C2 *= prime(prime(psidag(j), "Site"), jl);
           
            C3 *= psi(j);
            C3 *= op(sites, "LS+", j);
            C3 *= prime(prime(psidag(j), "Site"), jl);
            
            auto sij = elt(C1)+ elt(C2)+ elt(C3);
             
             SaveFile10<<i<<"     "<<j<<"     "<<sij<<endl;
             SaveFile10<<j<<"     "<<i<<"     "<<sij<<endl;
             
        }
    }
  SaveFile10.close();
  println("\n");



 //CdagC_nk
ofstream SaveFile11("nk.dat");
printfln("nk= ");    

 for (int i = 1; i <= N; i+=2)
   {     
        
   for (int j = i + 2; j <= N; j+=2)
      {
        psi.position(i);
        auto psidag = dag(psi);                 
        auto ir = commonIndex(psi(i), psi(i + 1), "Link");
        auto C1 = psi(i);
        C1 = noPrime(C1 * op(sites, "Adagup", i)) * op(sites, "F", i);
        C1 *=  prime(prime(psidag(i), "Site"), ir);

        auto C2 = psi(i);
        C2 *= op(sites, "Adagdn", i);
        C2 *= prime(prime(psidag(i), "Site"), ir);               

            for (int k = i + 1; k < j; ++k)
            {
                if(k%2 == 0) 
                {
                 C1 *= psi(k);
                 C1 *= prime(psidag(k), "Link"); 
                 C2 *= psi(k);
                 C2 *= prime(psidag(k), "Link"); 

                 }
               else
                {
              
                 C1 *= psi(k);  
                 C1 *= op(sites, "F", k); 
                 C1 *= prime(psidag(k));
                 C2 *= psi(k);  
                 C2 *= op(sites, "F", k); 
                 C2 *= prime(psidag(k));
                }
             }
            auto jl = commonIndex(psi(j), psi(j - 1), "Link");

            C1 *= psi(j);
            C1 *= op(sites, "Aup", j);
            C1 *= prime(prime(psidag(j), "Site"), jl);

            C2 *= psi(j);
            C2 = noPrime(C2 * op(sites, "F", j), "Site") * op(sites, "Adn", j);
            C2 *= prime(prime(psidag(j), "Site"), jl);

            auto cij1 = elt(C1); auto cij2 = elt(C2);
            auto cij = cij1 + cij2 ;
            SaveFile11<<i<<"     "<<j<<"     "<<cij<<endl;
            SaveFile11<<j<<"     "<<i<<"     "<<cij<<endl;
        
        }
    }
SaveFile11.close();
println("\n");

 //Nk
ofstream SaveFile12("Nk.dat");
printfln("Nk= ");    
 
 for (int i = 1; i <= N; i+=2)
   {     
        
   for (int j = i + 2; j <= N; j+=2)
      {
        psi.position(i);                
        auto ir = commonIndex(psi(i), psi(i + 1), "Link");
        auto C1 = psi(i);
        C1 *=  op(sites, "Ntot", i);
        C1 *=  dag(prime(prime(psi(i), "Site"), ir));

            for (int k = i + 1; k < j; ++k)
            {
                
                 C1 *= psi(k);
                 C1 *= dag(prime(psi(k), "Link"));    
               
            }
            auto jl = commonIndex(psi(j), psi(j - 1), "Link");
            C1 *= psi(j);
            C1 *= op(sites, "Ntot", j);
            C1 *= dag(prime(prime(psi(j), "Site"), jl));


        auto C2 = elt(psi(i)*op(sites,"Ntot",i)*dag(prime(psi(i),"Site")));
        
        psi.position(j);   
        auto C3 = elt(psi(j)*op(sites,"Ntot",j)*dag(prime(psi(j),"Site")));

        auto cij1 = C2*C3;
        auto cij = elt(C1)-cij1; 

        SaveFile12<<i<<"     "<<j<<"     "<<C2<<"     "<<C3<<"     "<<cij1<<"     "<<cij<<endl;
        SaveFile12<<j<<"     "<<i<<"     "<<C3<<"     "<<C2<<"     "<<cij1<<"     "<<cij<<endl;
          
    }
   }
SaveFile12.close();
println("\n");

*/


//Nk2
ofstream SaveFile13("Nk2.dat");
printfln("Nk2= ");    
 double n1[N+1],nn[N+1];
 printfln("iti_phixx= ");
   for (int i=0;i<N;i++) 
    {
      n1[i]=0;
	    nn[i]=0;
	
     }
 for (int i =N/4- 1; i <= 3*N/4; i+=2)
   {     
        
   for (int j = i + 2; j <= 3*N/4; j+=2)
      {
        psi.position(i);                
        auto ir = commonIndex(psi(i), psi(i + 1), "Link");
        auto C1 = psi(i);
        C1 *=  op(sites, "Ntot", i);
        C1 *=  dag(prime(prime(psi(i), "Site"), ir));

            for (int k = i + 1; k < j; ++k)
            {
                
                 C1 *= psi(k);
                 C1 *= dag(prime(psi(k), "Link"));    
               
            }
            auto jl = commonIndex(psi(j), psi(j - 1), "Link");
            C1 *= psi(j);
            C1 *= op(sites, "Ntot", j);
            C1 *= dag(prime(prime(psi(j), "Site"), jl));


        auto C2 = elt(psi(i)*op(sites,"Ntot",i)*dag(prime(psi(i),"Site")));
        
        psi.position(j);   
        auto C3 = elt(psi(j)*op(sites,"Ntot",j)*dag(prime(psi(j),"Site")));

        auto cij1 = C2*C3;
        auto cij = elt(C1)-cij1; 

           

        SaveFile13<<i<<"     "<<j<<"       "<<cij1<<"     "<<cij<<endl;
        SaveFile13<<j<<"     "<<i<<"       "<<cij1<<"     "<<cij<<endl;

        n1[(j-i)/2]++;
	      nn[(j-i)/2]+=cij;  
          
    }
   }

   for (int i=0;i<N;i++)
    {
      if (n1[i]!=0)
      SaveFile13<<i<<"      "<<nn[i]/n1[i]<<endl;
    }

SaveFile13.close();
println("\n");

/*
//pair correlation,iti  
 ofstream SaveFile6("iti_phixxr=1.dat");
 double corr1[N+1],corrs[N+1],corrt[N+1];
 double pairs,pairt;
 printfln("iti_phixx= ");
   for (int i=0;i<N;i++) 
    {
      corr1[i]=0;
	corrs[i]=0;
	corrt[i]=0;
     }

for (int i=81;i<N-79;i=i+2)
  {  
    
    psi.position(i);
    auto psidag=dag(psi);
    auto ir = commonIndex(psi(i), psi(i+1), "Link");
    auto C1 = psi(i);
    C1 = noPrime(C1*op(sites,"Adagup",i),"Site")*op(sites,"F",i);
    C1 *= prime(prime(psidag(i), "Site"), ir);  
    
        C1 *= psi(i+1);  
        C1 *= prime(psidag(i+1),"Link");
     
    C1 *= psi(i+2);
    C1 = noPrime(C1*op(sites,"Adagdn",i+2),"Site");
    C1 = noPrime(C1*op(sites,"F",i+2),"Site");

    //C1 *= prime(psidag(i+2));

    auto C2 = psi(i);
    C2 *= op(sites,"Adagdn",i);
    C2 *= prime(prime(psidag(i), "Site"), ir);
            
            C2 *= psi(i+1);  
            C2 *= prime(psidag(i+1),"Link");

    C2 *= psi(i+2);
    C2 = noPrime(C2*op(sites,"Adagup",i+2), "Site");
    //C2 *= prime(psidag(i+2));

      int b=2; int j=i+2;
      auto il=commonIndex(psi(j+b-1),psi(j+b),"Link");  
      auto D1 = noPrime(C1*op(sites,"F",j),"Site");
      D1 = noPrime(D1*op(sites,"Adn",j),"Site"); 
      D1 *= op(sites,"F",j);
      D1 *= prime(psidag(j));

            D1 *= psi(j+1);  
            D1 *= prime(psidag(j+1),"Link");
           
      D1 *= psi(j+b);
      D1 *= op(sites,"Aup",j+b);
      D1 *= prime(prime(psidag(j+b),"Site"),il);

      auto D2 = noPrime(C1*op(sites,"Aup",j),"Site");
      D2 *= op(sites,"F",j);
	D2 *= prime(psidag(j));
      
            D2 *= psi(j+1);  
             D2 *= prime(psidag(j+1),"Link");	

      D2 = noPrime(D2*psi(j+b)*op(sites,"F",j+b),"Site");          
      D2 *= op(sites,"Adn",j+b);
	D2 *= prime(prime(psidag(j+b),"Site"),il);
  
      pairs = elt(D1) - elt(D2);
      pairt = elt(D1) + elt(D2);
      
      auto D3 = noPrime(C2*op(sites,"F",j),"Site");
      D3 = noPrime(D3*op(sites,"Adn",j),"Site"); 
      D3 = D3*op(sites,"F",j);
      D3 = D3*prime(psidag(j));
      
            D3 *= psi(j+1);  
            D3 *= prime(psidag(j+1),"Link");
                     
      D3 *= psi(j+b);
      D3 *= op(sites,"Aup",j+b);
      D3 *= prime(prime(psidag(j+b),"Site"),il);

	auto D4 = noPrime(C2*op(sites,"Aup",j),"Site");
      D4 *= op(sites,"F",j);
	D4 *= prime(psidag(j));
	
            D4 *= psi(j+1);  
            D4 *= prime(psidag(j+1),"Link");

	D4 = noPrime(D4*psi(j+b)*op(sites,"F",j+b),"Site");
      D4 *= op(sites,"Adn",j+b);
	D4 *= prime(prime(psidag(j+b),"Site"),il);

      pairs = pairs - elt(D3) + elt(D4);
	pairt = pairt + elt(D3) + elt(D4);

	//printf("%f  %f",elt(D3),elt(D4));
	//printf("  %d  %d  %.10f  %.10f\n",i,j,pairs,pairt);

      SaveFile6<<elt(D1)<<"  "<<elt(D2)<<"  "<<elt(D3)<<"  "<<elt(D4)<<"  "<<i<<"  "<<j<<"  "<<pairs<<"  "<<pairt<<endl;
       
      corr1[(j-i)/2]++;
	corrs[(j-i)/2]+=pairs;
	corrt[(j-i)/2]+=pairt;
		
          
    }
          
 for (int i=0;i<N;i++)
  if (corr1[i]!=0)

SaveFile6<<i<<"     "<<corrs[i]/corr1[i]<<"     "<<corrt[i]/corr1[i]<<endl;
SaveFile6.close();
println("\n");  




//pair correlation,iti  
 ofstream SaveFile7("iti_phixx.dat");
  printfln("iti_phixx= ");
   for (int i=0;i<N;i++) 
    {
      corr1[i]=0;
	corrs[i]=0;
	corrt[i]=0;
     }

for (int i=81;i<N-79;i=i+2)
  {  
    int a=2; int b=2; 
    psi.position(i);
    auto psidag=dag(psi);
    auto ir = commonIndex(psi(i), psi(i+1), "Link");
    auto C1 = psi(i);
    C1 = noPrime(C1*op(sites,"Adagup",i),"Site")*op(sites,"F",i);
    C1 *= prime(prime(psidag(i), "Site"), ir);  
    for (int inter=i+1;inter<=i+a-1;inter++)
      { 
        C1 *= psi(inter);  
        C1 *= prime(psidag(inter),"Link");

       }      
    C1 *= psi(i+a);
    C1 = noPrime(C1*op(sites,"Adagdn",i+a),"Site")*op(sites,"F",i+a);   
    C1 *= prime(psidag(i+a));

    auto C2 = psi(i);
    C2 *= op(sites,"Adagdn",i);
    C2 *= prime(prime(psidag(i), "Site"), ir);
    for (int inter=i+1;inter<=i+a-1;inter++)
	  {
            C2 *= psi(inter);  
            C2 *= prime(psidag(inter),"Link");

        }
    C2 *= psi(i+a)*op(sites,"Adagup",i+a);
    C2 *= prime(psidag(i+a));
 
 for (int j=i+a+1;j<N-79;j++)
    {
   
   if (j%2 == 1) 
     { 
      auto il=commonIndex(psi(j+b-1),psi(j+b),"Link");  
      auto D1 = noPrime(C1*psi(j)*op(sites,"F",j),"Site");
      D1 = noPrime(D1*op(sites,"Adn",j),"Site"); 
      D1 *= op(sites,"F",j);
      D1 *= prime(psidag(j));
      for (int inter=j+1;inter<=j+b-1;inter++) 
         {
            D1 *= psi(inter);  
            D1 *= prime(psidag(inter),"Link");

          }           
      D1 *= psi(j+b);
      D1 *= op(sites,"Aup",j+b);
      D1 *= prime(prime(psidag(j+b),"Site"),il);

      auto D2 = noPrime(C1*psi(j)*op(sites,"Aup",j),"Site");
      D2 *= op(sites,"F",j);
	D2 *= prime(psidag(j));
      for (int inter=j+1;inter<=j+b-1;inter++)
	   {
            D2 *= psi(inter);  
            D2 *= prime(psidag(inter),"Link");

         }
	D2 = noPrime(D2*psi(j+b)*op(sites,"F",j+b),"Site");          
      D2 *= op(sites,"Adn",j+b);
	D2 *= prime(prime(psidag(j+b),"Site"),il);
  
      pairs = elt(D1) - elt(D2);
      pairt = elt(D1) + elt(D2);
      
      auto D3 = noPrime(C2*psi(j)*op(sites,"F",j),"Site");
      D3 = noPrime(D3*op(sites,"Adn",j),"Site"); 
      D3 = D3*op(sites,"F",j);
      D3 = D3*prime(psidag(j));
      for (int inter=j+1;inter<=j+b-1;inter++) 
         {
            D3 *= psi(inter);  
            D3 *= prime(psidag(inter),"Link");
          }           
      D3 *= psi(j+b);
      D3 *= op(sites,"Aup",j+b);
      D3 *= prime(prime(psidag(j+b),"Site"),il);

	auto D4 = noPrime(C2*psi(j)*op(sites,"Aup",j),"Site");
      D4 *= op(sites,"F",j);
	D4 *= prime(psidag(j));
	for (int inter=j+1;inter<=j+b-1;inter++)
	   {
            D4 *= psi(inter);  
            D4 *= prime(psidag(inter),"Link");

          }
	D4 = noPrime(D4*psi(j+b)*op(sites,"F",j+b),"Site");
      D4 *= op(sites,"Adn",j+b);
	D4 *= prime(prime(psidag(j+b),"Site"),il);

      pairs = pairs - elt(D3) + elt(D4);
	pairt = pairt + elt(D3) + elt(D4);

	//printf("%f  %f",elt(D3),elt(D4));
	//printf("  %d  %d  %.10f  %.10f\n",i,j,pairs,pairt);

      SaveFile7<<elt(D1)<<"  "<<elt(D2)<<"  "<<elt(D3)<<"  "<<elt(D4)<<"  "<<i<<"  "<<j<<"  "<<pairs<<"  "<<pairt<<endl;
       
  corr1[(j-i)/2]++;
	corrs[(j-i)/2]+=pairs;
	corrt[(j-i)/2]+=pairt;
		
        }
          C1 *= psi(j);
          C1 *= prime(psidag(j),"Link");
          C2 *= psi(j);
          C2 *= prime(psidag(j),"Link");     

      }
   }       
 for (int i=0;i<N;i++)
  if (corr1[i]!=0)

SaveFile7<<i<<"     "<<corrs[i]/corr1[i]<<"     "<<corrt[i]/corr1[i]<<endl;
SaveFile7.close();
println("\n");  



//Ntot
ofstream SaveFile14("iti_Ntot.dat");
printfln("iti_Ntot= ");
  for (int i=1;i<N;i=i+2)
     {
       psi.position(i);
       auto psidag = dag(psi);                 
       auto C = elt(psi(i)*op(sites,"Ntot",i)*prime(psidag(i),"Site"));
       SaveFile14<<i<<"     "<<C<<endl; 
       
      }
  SaveFile14.close();
  println("\n");


//dd
ofstream SaveFile15("iti_dd.dat");
printfln("iti_dd= ");
  for (int i=1;i<N;i=i+2)
     {
       psi.position(i);
       auto psidag = dag(psi);                 
       auto C = elt(psi(i)*op(sites,"Nupdn",i)*prime(psidag(i),"Site"));
       SaveFile15<<i<<"     "<<C<<endl; 
       
      }
  SaveFile15.close();
  println("\n");



//nj,itinerant
ofstream SaveFile18("iti_nj.dat");
printfln("iti_nj = ");
for (int i = 1; i < N; i+=2)
   {
        psi.position(i);
        auto psidag = dag(psi);                 
        auto C1 = psi(i);
        C1 = noPrime(C1 * op(sites, "Aup", i)) * op(sites, "Adagup", i);
        C1 *=  prime(psidag(i), "Site");

        auto C2 = psi(i);
        C2 = noPrime(C2 * op(sites, "Adn", i)) * op(sites, "Adagdn", i);
        C2 *= prime(psidag(i), "Site");            

         auto cij1 = elt(C1); auto cij2 = elt(C2);
         auto cij = cij1 + cij2 ;
            
          SaveFile18<<i<<"    "<<cij1<<"    "<<cij2<<"    "<<cij<<endl;
          
       }
    
    SaveFile18.close();
    println("\n");  
*/


  return 0;
}