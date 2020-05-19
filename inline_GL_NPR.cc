#include "inline_Ob.h"
#include "chroma.h"
#include <math.h>
#include <complex>
#include <iostream>
#include <string>
#include "stdio.h"



namespace Chroma 
{
    
    namespace InlineObEnv 
    {
	//Name of the measurement to be called in the XML input file
	const std::string name = "GMF_O_b";
	
	//This function is used with the factory thing
	AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
						const std::string& path) 
	{
	    //Create new instance of measurement class using params
	    //from the passed xml file
	    return new InlineMyMeas(InlineObParams(xml_in, path));
	}
		
	// Local registration flag
	namespace {
	    bool registered = false;
	}
	
	// Function to register all the factories
	bool registerAll() 
	{
	    bool success = true; 
	    if (! registered)
	    {
		success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
		QDPIO::cout << "Registering " << name << " " << success << std::endl;
		registered = true;
	    }
	    return success;
	}
	
	
	
	/*** Implementation of Parameter functions ***/
	
	//Set default parameters
	InlineObParams::InlineObParams() { frequency = 0; radius = 0; srcs.resize(1); }
	
	//Read parameters in from xml file
	InlineObParams::InlineObParams(XMLReader& xml_in, const std::string& path) 
	{
	    try 
	    {
		XMLReader paramtop(xml_in, path);
		
		if (paramtop.count("Frequency") == 1)
		    read(paramtop, "Frequency", frequency);
		else
		    frequency = 1;
		
		read(paramtop, "NamedObject", named_obj);

		//Read in the starting position time (ie the time
		//of the first source) and the time between each
		//source. This assumes that each source is equally
		//spaced in time. See branch First-O_b for code
		//that drops this assumtion and takes sources
		//as the arguments (loc & start/stop times)
		read(paramtop, "Multi_Src", srcs);

		if(paramtop.count("radius") == 1)
		    read(paramtop, "radius", radius);
		else
		    radius = 0;

		
		// Possible alternate XML file pattern
		if (paramtop.count("xml_file") != 0) 
		{
		    read(paramtop, "xml_file", xml_file);
		} else
		    xml_file = "";
		
	    }
	    catch(const std::string& e) 
	    {
		QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
		QDP_abort(1);
	    }
	}
	
	
	// Write loaded params to xml file
	void InlineObParams::write(XMLWriter& xml_out, const std::string& path) 
	{
	    push(xml_out, path);
	    
	    // write our all params
	    QDP::write(xml_out, "Multi_Src", srcs);
	    //QDP::write(xml_out, "Named_Object", named_obj);
	    QDP::write(xml_out, "radius", radius);

		  
	    
	    if(xml_file != "")
		QDP::write(xml_out, "xml_file", xml_file);
	    pop(xml_out);
	}
    } //End namespace InlineObEnv

    /*** Inline Measurement function implimentation ***/
    // Function call
    void InlineMyMeas::operator()(unsigned long update_no,
				  XMLWriter& xml_out) 
    {
	// If xml file not empty, then use alternate
	if (params.xml_file != "")
	{
	    std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	    push(xml_out, "GMF_O_b");
	    write(xml_out, "update_no", update_no);
	    write(xml_out, "xml_file", xml_file);
	    pop(xml_out);
	    
	    XMLFileWriter xml(xml_file);
	    func(update_no, xml);
	}
	else
	{
	    func(update_no, xml_out);
	}
    }
    
    /* This function is used to obtain the gluon field F_\mu\nu shift "len" times along "dir" direction.
     * The input F is one matrix element of F_\mu\nu.
     */
    LatticeColorMatrix  field(int dir, int len, LatticeColorMatrix F)
    {
        LatticeColorMatrix Ftem, Fshift; /* Ftem is the temporary variable,
                                            Fshift is the gluon field F_\mu\nu shift "len" times along "dir" direction that will be returned.
                                          */
        Fshift = F;
        for(int i = 1; i <= len; i++)
        {
                Ftem = Fshift;
                Fshift = shift(Ftem, FORWARD, dir);
        }
        return Fshift;
    }
 
    /* This function is used to construct and return the wilson line at "dir" direction with the length "len".
     * The input u is the unit wilson link at "dir" direction.
     * U(x,x+len\hat{dir})= \sum_{y=x}^{x+len-1}U(y,y+1\hat{dir})
     */
    LatticeColorMatrix  wilsonline(int dir, int len, LatticeColorMatrix u) 
    {
	LatticeColorMatrix utem, ushift, uline; /* utem is the temporary wilson line
						ushift is the shift "len" times wilson link
						uline is the wilson line at "dir" direction with the length "len" that will be returned
						*/
	ushift = u;
	uline = u;
	for(int i = 1; i < len; i++)
        {
		utem = ushift;
		ushift = shift(utem, FORWARD, dir);
		utem = uline;
		uline = utem * ushift;
	}
	return uline; 
    }

    /* The function is to calculte the m*n plaquette at the direction of
     * \mu and \nu. The unit wilson link umu and unu at \mu and \nu direction are
     * necessary input.
     * P1_{\mu\nu}=U(x,x+m\hat{\mu})U(x+m\hat{\mu},x+m\hat{\mu}+n\hat{\nu})U(x+m\hat{\mu}+n\hat{\nu},x+n\hat{\nu})U(x+n\hat{\nu},x)
     * P2_{\mu\nu}=P1_{\nu,-\mu}
     * P3_{\mu\nu}=P1_{-\mu,-\nu}
     * P4_{\mu\nu}=P1_{-\nu,\mu}
     */ 
    LatticeColorMatrix  plaquette(int mu, int nu, int m, int n, LatticeColorMatrix umu, LatticeColorMatrix unu, int pla_num)
    {
        LatticeColorMatrix wltem, wlmu, wlnu, wlmu_shift, wlnu_shift, plane_plaq_mn; /* wltem is the intermediated line
											wlmu is the first quater of the wilson line of the plaquette
                                                                                        wlnu is the second quater of the wilson line of the plaquette
                                                                                        wlmu_shift is the third quater of the wilson line of the plaquette
                                                                                        wlnu_shift is the last quater of the wilson line of the plaquette
											plane_plaq_mn is the m*n plaquette that will be returned
											*/
	wlmu = wilsonline(mu, m, umu); 
	wlnu = wilsonline(nu, n, unu);
	wlmu_shift = wlmu;
        wlnu_shift = wlnu;
        for(int i = 1; i <= n; i++)
        {
                wltem = wlmu_shift;
                wlmu_shift = shift(wltem, FORWARD, nu);
        }
        for(int j = 1; j <= m; j++)
        {
                wltem = wlnu_shift;
                wlnu_shift = shift(wltem, FORWARD, mu);
        }
	switch (pla_num)  //generate 4 plaquette used in the construction of F_\mu\nu
        {
                case 1:plane_plaq_mn = wlmu*wlnu_shift*adj(wlmu_shift)*adj(wlnu);
                break;
                case 2:plane_plaq_mn = shift(wlnu_shift*adj(wlmu_shift)*adj(wlnu)*wlmu, BACKWARD, mu);
                break;
                case 3:plane_plaq_mn = shift(shift(adj(wlmu_shift)*adj(wlnu)*wlmu*wlnu_shift, BACKWARD, mu), BACKWARD, nu);
                break;
                case 4:plane_plaq_mn = shift(adj(wlnu)*wlmu*wlnu_shift*adj(wlmu_shift), BACKWARD, nu);
                break;
        }
	//plane_plaq_mn = ulinem*ulinenshift*adj(ulinemshift)*adj(ulinen);
	return plane_plaq_mn;
    }

    /* The function is used to calculate the mean(real(trace())) of the plaquette.
     * The input should be a matrix element of plane plaquette.
     * The output is a matrix element of the mean(real(trace())) of the plaquette.
     */ 
    Double tr_plane_pla(LatticeColorMatrix plane_plaq)
    {
	Double tr_plane_plaq;
	tr_plane_plaq = sum(real(trace(plane_plaq)));
	tr_plane_plaq /= Double(Layout::vol() * Nc); //average on sites and colors
	tr_plane_plaq = tr_plane_plaq; //symmetric
	return tr_plane_plaq;
    }


    /* The function is used to construct the O_3 operator (both local and non-local operator) and its 16 components.
     * The inputs are the wilson line length "len" and the direction "dir", 4 gluon field components (plaq[i]-adj(plaq[i])), and the unit wilson link u
     * The outputs are the O_3 16 componets at 3 direction (x,y,t if dir=z) and O_3 itself.
     */

    int delta(int i, int j)
    {
	if(i==j)
		delta=1;
	else
		delta=0;
	return delta;
    }

    Double lambda_tree(int rep, int alpha, int beta, int tau, multi1d<Double>)
    {
	if(rep==3)
		lambda_tree=(Nc*Nc-1)/sqrt(2)*(p[alpha]*p[alpha]+p[beta]*p[beta]);
	else
		lambda_tree=(Nc*Nc-1)/sqrt(2)*(2*(-I)*delta(beta,4)*p[alpha]*p[alpha]);
	return lambda_tree;	
    }

    

    Double Z_RIMOM(int alpha, int beta, int tau, LatticeColorMatrix O, LatticeColorMatrix u,multi2d<Double> p, multi2d<Double> x)
    {
	LatticeColorMatrix A_x, A_p, A_mp;
	Doulbe G_3pt, G_2pt;
	A_x=1/(2*I*g0)*((u-adj(u)-1/Nc*trace(u-adj(u))));
	for(i=0; i<Nc; i++)
		p_dot_x=p_dot_x+p[i]*x[i];
	A_p=sumx(exp(-I*p_dot_x)*A_x);
	A_mp=sumx(exp(I*p_dot_x)*A_x);
	G_3pt=sum(trace(O)*trace(A_p*A_mp))/(Layout::vol());
	G_2pt=sum(trace(A_p*A_mp))/(Layout::vol());
	Z_RIMOM=4*p2*G_3pt/lambda_tree[rep,nu,alpha,beta,tau,p]/G_2pt
	return Z_RIMOM;
    }


    multi2d<LatticeColorMatrix> fun_Opcomp(int len, int dir, multi2d<LatticeColorMatrix> F, LatticeColorMatrix u)
    {
        LatticeColorMatrix un;
        multi2d<LatticeColorMatrix> Fn;
        Fn.resize(Nd,Nd);
        Fn = 0;
        multi2d<LatticeColorMatrix> Opcomp;
        Opcomp.resize(3,Nd);
        Opcomp = 0;
        un = wilsonline(dir, len, u); //calculate the wilson line

        for(int mu = 0; mu < Nd; mu++)
        {
                for(int nu = mu+1; nu < Nd; nu++)
                {
                        Fn[nu][mu] = field(dir, len, F[nu][mu]); //shift the F_\mu\nu
                        Fn[mu][nu] = -Fn[nu][mu];  //anti-symmetric
                }
        }


        if(len==0) //local operators
        {
                for(int i=0;i<4; i++)
                {
                        Opcomp[0][i] = F[dir][i]*F[dir][i];
			Opcomp[1][i] = F[3][i]*F[dir][i];
			Opcomp[2][i] = F[3][i]*F[3][i];	       
                }
        }


        else //non-local
        {
                for(int i=0;i<4; i++)
                {
                        Opcomp[0][i] = F[dir][i]*un*F[dir][i]*adj(un);
                        Opcomp[1][i] = F[3][i]*un*F[dir][i]*adj(un);
                        Opcomp[2][i] = F[3][i]*un*F[3][i]*adj(un);
                }
        }
	
	return Opcomp;
    }





    /* The function is used to construct the operators (both local and non-local operators). 
     * The inputs are the wilson line length "len" and the direction "dir", gluon field F_\mu\nu, and the unit wilson link u
     * The outputs are the 6 operators.
     */ 
    multi1d<LatticeColorMatrix> fun_Operator(int len, int dir, multi2d<LatticeColorMatrix> F, LatticeColorMatrix u)
    {
	LatticeColorMatrix un;
	multi2d<LatticeColorMatrix> Fn;
	Fn.resize(Nd,Nd);
	Fn = 0;
	multi1d<LatticeColorMatrix> Op;
	Op.resize(11);
	Op = 0;

	un = wilsonline(dir, len, u); //calculate the wilson line

	for(int mu = 0; mu < Nd; mu++)
        {
        	for(int nu = mu+1; nu < Nd; nu++)
                {
                	Fn[nu][mu] = field(dir, len, F[nu][mu]); //shift the F_\mu\nu
                        //Fn[mu][nu] = Fn[nu][mu];  //symmetric
			Fn[mu][nu] = -Fn[nu][mu];  //anti-symmetric
                }
        }

/*
        LatticeColorMatrix um;
	un = adj(un);
        multi2d<LatticeColorMatrix> Fm;
        Fm.resize(Nd,Nd);
	for(int k=0;k<n; k++)
	{
		um = un;
		un = shift(um, BACKWARD, z);
	}
	for(int mu = 0; mu < Nd; mu++)
        {
                for(int nu = mu+1; nu < Nd; nu++)
                {
			Fm[nu][mu] = F[nu][mu];
                        for(int k=0;k<n; k++)
			{
				Fn[nu][mu] = shift(Fm[nu][mu], BACKWARD, z);
				Fm[nu][mu] = Fn[nu][mu];
			}
                        Fn[mu][nu] = -Fn[nu][mu];  //anti-symmetric
                }
        }
*/	



	/* Operators definition
 	 * O(F_{\mu \nu},F^{\alpha \beta},len) = tr(F_{\mu \nu}(0)U(0,len)F_{\mu \nu}(len)U(0,len)), \mu = 0,1,2,3, i=0,1,2 (without dir)
 	 * ZY Fan at el. Physical review letters, 121(24), 242001.
 	 * Op[0]: O_0 = O(F_{\mu t},F^{\mu t},len)-1/4 O(F_{\mu \nu},F^{\mu \mu},len)
 	 * Op[1]: O_1 = O(F_{\mu t},F^{\mu dir},len)
 	 * Op[2]: O_2 = O(F_{\mu dir},F^{\mu dir},len)-1/4 O(F_{\mu \nu},F^{\mu \mu},len)
 	 * Op[3]: O_3 = O(F_{\mu dir},F^{\mu dir},len)
 	 * Jianhui Zhang at el. 10.1103/PhysRevLett.122.142001
 	 * Op[6]: O_1 = O(F^{t i},F^{t i},len)
 	 * Op[7]: O_2 = O(F^{dir i},F^{dir i},len)
 	 * Op[5]: O_3 = O(F^{t i},F^{dir i},len)
 	 * Op[3]: O_4 = O(F^{dir \mu},F^{dir \mu},len)
 	 * Op[9]: \delta O_1 = \epsilon_ij O(F^{t i},F^{t j},len)
 	 * Op[8]: \delta O_2 = \epsilon_ij O(F^{dir i},F^{dir j},len)
 	 * Op[10]: \delta O_3 = \epsilon_ij O(F^{t i},F^{dir j},len)
 	 * arXiv:1910.13963v1
 	 * Op[4]: O = O(F^{dir i},F^{t i},len)
 	 */

/*      //This is a test on the wrong combination of the links 
        LatticeColorMatrix um, ushift;
        un = u;
	ushift = u;
        multi2d<LatticeColorMatrix> Fm;
        Fm.resize(Nd,Nd);
        for(int k=1;k<len; k++)
        {
                um = ushift;
		ushift = shift(um, FORWARD, dir);
		um = un;
                un = ushift*um;
        }
*/	
	if(len==0) //local operators
	{
		for(int i=0;i<4; i++)
                { 
			Op[0] += F[3][i]*F[3][i];
			Op[1] += F[3][i]*F[dir][i];
			Op[2] += F[dir][i]*F[dir][i];
                	Op[3] += F[dir][i]*F[dir][i];
			//Op[4] += F[3][i]*F[3][i];
			
			//Op[5] += F[dir][i]*F[dir][i];
                }
		for(int i=0;i<4; i++)
                {
                        for(int j=0;j<i;j++)
                        {
                                Op[0] -= 0.5*F[j][i]*F[j][i];
                                Op[2] -= 0.5*F[j][i]*F[j][i];
                        }
               	}
                for(int i=0;i<4; i++)
                {
			if(i!=dir && i!=3)
			{
				Op[4] += F[dir][i]*F[3][i];
				Op[5] += F[3][i]*F[dir][i];
                                Op[6] += F[3][i]*F[3][i];
                                Op[7] += F[dir][i]*F[dir][i];
			}
		//Op[6] = F[dir][0]*F[dir][0];
		//Op[7] = F[dir][1]*F[dir][1];
		//Op[8] = F[dir][2]*F[dir][2];
		//Op[9] = F[dir][3]*F[dir][3];
		}
                if(dir==0)
                {
                        Op[8] = F[3][1]*F[3][2]-F[3][2]*F[3][1];
                        Op[9] = F[dir][1]*F[dir][2]-F[dir][2]*F[dir][1];
                        Op[10] = F[3][1]*F[dir][2]-F[3][2]*F[dir][1];
                }
                else if(dir==1)
                {
                        Op[8] = F[3][2]*F[3][0]-F[3][0]*F[3][2];
                        Op[9] = F[dir][2]*F[dir][0]-F[dir][0]*F[dir][2];
                        Op[10] = F[3][2]*F[dir][0]-F[3][0]*F[dir][2];
                }
                else if(dir==2)
                {
                        Op[8] = F[3][0]*F[3][1]-F[3][1]*F[3][0];
                        Op[9] = F[dir][0]*F[dir][1]-F[dir][dir]*F[3][0];
                        Op[10] = F[3][0]*F[dir][1]-F[3][1]*F[dir][0];
                }

 
	}

	else //non-local
	{

		for(int i=0;i<4; i++)
                {
			
			
                        Op[0] += F[3][i]*un*Fn[3][i]*adj(un);
                        Op[1] += F[3][i]*un*Fn[dir][i]*adj(un);
                        Op[2] += F[dir][i]*un*Fn[dir][i]*adj(un);
                        Op[3] += F[dir][i]*un*Fn[dir][i]*adj(un);
                        //Op[4] += F[3][i]*un*Fn[3][i]*adj(un);
			
                        
                        //Op[5] += F[dir][i]*un*Fn[dir][i]*adj(un);
                }
                for(int i=0;i<4; i++)
                {
                        for(int j=0;j<i;j++)
                        {
                                Op[0] -= 0.5*F[j][i]*un*Fn[j][i]*adj(un);
                                Op[2] -= 0.5*F[j][i]*un*Fn[j][i]*adj(un);
                        }
                }
                for(int i=0;i<4; i++)
                {
                        if(i!=dir && i!=3)
                        {
				Op[4] += F[dir][i]*un*Fn[3][i]*adj(un);
				Op[5] += F[3][i]*un*Fn[dir][i]*adj(un);
                                Op[6] += F[3][i]*un*Fn[3][i]*adj(un);
                                Op[7] += F[dir][i]*un*Fn[dir][i]*adj(un);
                        }
		}
		
		if(dir==0)
		{	
			Op[8] = F[3][1]*un*Fn[3][2]*adj(un)-F[3][2]*un*Fn[3][1]*adj(un);
			Op[9] = F[dir][1]*un*Fn[dir][2]*adj(un)-F[dir][2]*un*Fn[dir][1]*adj(un);
			Op[10] = F[3][1]*un*Fn[dir][2]*adj(un)-F[3][2]*un*Fn[dir][1]*adj(un);
		}
		else if(dir==1)
		{
                        Op[8] = F[3][2]*un*Fn[3][0]*adj(un)-F[3][0]*un*Fn[3][2]*adj(un);
                        Op[9] = F[dir][2]*un*Fn[dir][0]*adj(un)-F[dir][0]*un*Fn[dir][2]*adj(un);
                        Op[10] = F[3][2]*un*Fn[dir][0]*adj(un)-F[3][0]*un*Fn[dir][2]*adj(un);
		}
		else if(dir==2)
                {
                        Op[8] = F[3][0]*un*Fn[3][1]*adj(un)-F[3][1]*un*Fn[3][0]*adj(un);
                        Op[9] = F[dir][0]*un*Fn[dir][1]*adj(un)-F[dir][1]*un*Fn[dir][0]*adj(un);
                        Op[10] = F[3][0]*un*Fn[dir][1]*adj(un)-F[3][1]*un*Fn[dir][0]*adj(un);
                }


		//Op[6] += F[dir][0]*un*Fn[dir][0]*adj(un);
                //Op[7] += F[dir][1]*un*Fn[dir][1]*adj(un);
                //Op[8] += F[dir][2]*un*Fn[dir][2]*adj(un);
                //Op[9] += F[dir][3]*un*Fn[dir][3]*adj(un);		
/*
                for(int i=0;i<4; i++)
                {
                        Op[0] += Fn[3][i]*adj(un)*F[3][i]*un;
                        Op[1] += Fn[3][i]*adj(un)*F[z][i]*un;
                        Op[2] += Fn[z][i]*adj(un)*F[z][i]*un;
                        Op[3] += Fn[z][i]*adj(un)*F[z][i]*un;
                        Op[4] += Fn[3][i]*adj(un)*F[3][i]*un;
                        if(i==z || i==3)
                        { Op[5] += Fn[z][i]*adj(un)*F[z][i]*un;}
                }
                for(int i=0;i<4; i++)
                {
                        for(int j=0;j<i;j++)
                        {
                                Op[0] -= 0.5*Fn[j][i]*adj(un)*F[j][i]*un;
                                Op[2] -= 0.5*Fn[j][i]*adj(un)*F[j][i]*un;
                        }
                }
*/


	}

	return Op;
	
    }



    /*** Measurement code stars here ***/
    void InlineMyMeas::func(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
	START_CODE();
	
	QDPIO::cout << InlineObEnv::name << ": Begining" << std::endl;
	
	StopWatch snoop;
	snoop.reset();
	snoop.start();

	//Print out boilerplate stuff to xml
	push(xml_out, "GMF_O_b");
	write(xml_out, "update_no", update_no);

	//Write out the input
	params.write(xml_out, "Input");


	/** Calculate the two dimensional plaquettes **/
	
	//Get link matrix
	multi1d<LatticeColorMatrix> u;
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);	
	//chroma function to calculate the plaquette
	Wloop(xml_out, "GMF_O_b", u);
	//Create variables to store the plaquette and the average
	//trace of the plaquette
        multi3d<LatticeColorMatrix> plane_plaq;
        plane_plaq.resize(4,Nd,Nd); //Quadrent, plane

        /* Calculate the plaquette and the average trace of
         * the plaquette
         */
        for(int mu = 0; mu < Nd; mu++)
        {
            for(int nu = mu+1; nu < Nd; nu++)
            {
                plane_plaq[0][nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 1);
		plane_plaq[1][nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 2);
                plane_plaq[2][nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 3);
                plane_plaq[3][nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 4);
		//plane_plaq[0][mu][nu] = plane_plaq[0][nu][mu];
	    }
	}

        
	multi2d<LatticeColorMatrix> plane_plaq_11, plane_plaq_12, plane_plaq_13, plane_plaq_14, plane_plaq_22, plane_plaq_32; //plane_plaq_{m}{n} represents the n*m plaquette
        multi2d<Double> tr_plane_plaq_11, tr_plane_plaq_12, tr_plane_plaq_13, tr_plane_plaq_14, tr_plane_plaq_22, tr_plane_plaq_32; //The mean(trace(real of the plaquette
        plane_plaq_11.resize(Nd,Nd);
        plane_plaq_12.resize(Nd,Nd);
        plane_plaq_13.resize(Nd,Nd);
        plane_plaq_14.resize(Nd,Nd);
        plane_plaq_22.resize(Nd,Nd);
        plane_plaq_32.resize(Nd,Nd);
        tr_plane_plaq_11.resize(Nd,Nd);
        tr_plane_plaq_12.resize(Nd,Nd);
        tr_plane_plaq_13.resize(Nd,Nd);
        tr_plane_plaq_14.resize(Nd,Nd);
        tr_plane_plaq_22.resize(Nd,Nd);
        tr_plane_plaq_32.resize(Nd,Nd);

        for(int mu = 0; mu < Nd; mu++)
        {
        	for(int nu = mu+1; nu < Nd; nu++)
                {
			plane_plaq_11[nu][mu] = plaquette(mu, nu, 1, 1, u[mu], u[nu], 1);
			plane_plaq_12[nu][mu] = plaquette(mu, nu, 1, 2, u[mu], u[nu], 1);
                       	plane_plaq_13[nu][mu] = plaquette(mu, nu, 1, 3, u[mu], u[nu], 1);
                        plane_plaq_14[nu][mu] = plaquette(mu, nu, 1, 4, u[mu], u[nu], 1);
                       	plane_plaq_22[nu][mu] = plaquette(mu, nu, 2, 2, u[mu], u[nu], 1);
                        plane_plaq_32[nu][mu] = plaquette(mu, nu, 3, 2, u[mu], u[nu], 1);
                }
        }
	for(int mu = 0; mu < Nd; mu++)
	{
		for(int nu = mu+1; nu < Nd; nu++)
		{
			tr_plane_plaq_11[mu][nu] = tr_plane_pla(plane_plaq_11[nu][mu]);
                        tr_plane_plaq_12[mu][nu] = tr_plane_pla(plane_plaq_12[nu][mu]);
                        tr_plane_plaq_13[mu][nu] = tr_plane_pla(plane_plaq_13[nu][mu]);
                        tr_plane_plaq_14[mu][nu] = tr_plane_pla(plane_plaq_14[nu][mu]);
                        tr_plane_plaq_22[mu][nu] = tr_plane_pla(plane_plaq_22[nu][mu]);
                        tr_plane_plaq_32[mu][nu] = tr_plane_pla(plane_plaq_32[nu][mu]);
		}
	}	
	
	//Pirnt the mean(trace(real of the plaquette.
        for(int mu = 0; mu < Nd; mu++)
        {
        	for(int nu = mu+1; nu < Nd; nu++)
        	{
                	write(xml_out, "plane_plaq_11_" + std::to_string(mu) +
                	std::to_string(nu), tr_plane_plaq_11[mu][nu]);
                }
        }

        for(int mu = 0; mu < Nd; mu++)
        {
                for(int nu = mu+1; nu < Nd; nu++)
                {
                        write(xml_out, "plane_plaq_12_" + std::to_string(mu) +
                        std::to_string(nu), tr_plane_plaq_12[mu][nu]);
                }
        }

        for(int mu = 0; mu < Nd; mu++)
        {
                for(int nu = mu+1; nu < Nd; nu++)
                {
                        write(xml_out, "plane_plaq_13_" + std::to_string(mu) +
                        std::to_string(nu), tr_plane_plaq_13[mu][nu]);
                }
        }

        for(int mu = 0; mu < Nd; mu++)
        {
                for(int nu = mu+1; nu < Nd; nu++)
                {
                        write(xml_out, "plane_plaq_14_" + std::to_string(mu) +
                        std::to_string(nu), tr_plane_plaq_14[mu][nu]);
                }
        }

        for(int mu = 0; mu < Nd; mu++)
        {
                for(int nu = mu+1; nu < Nd; nu++)
                {
                        write(xml_out, "plane_plaq_22_" + std::to_string(mu) +
                        std::to_string(nu), tr_plane_plaq_22[mu][nu]);
                }
        }

        for(int mu = 0; mu < Nd; mu++)
        {
                for(int nu = mu+1; nu < Nd; nu++)
                {
                        write(xml_out, "plane_plaq_32_" + std::to_string(mu) +
                        std::to_string(nu), tr_plane_plaq_32[mu][nu]);
                }
        }

	//loop over direction 0,1,2
	for(int dir = 0; dir < 3; dir++)
        {
                QDPIO::cout << "Finding F" << std::endl;
                int mn = 11; //maximum wilson length
                /** Find F_{n,munu} **/
                multi2d<LatticeColorMatrix> F, Fn;
                F.resize(Nd,Nd);
                Fn.resize(Nd,Nd);
		F = 0;
		Fn = 0;

                // Calculate the F
                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {
                        for(int i = 0; i < 4; i++)
                        {
                                F[nu][mu] += plane_plaq[i][nu][mu] - adj(plane_plaq[i][nu][mu]);
                        }
                        F[nu][mu] *= 1/8.0;
                        //F[mu][nu]= F[nu][mu]; //symmetric
			F[mu][nu]= -F[nu][mu]; //anti-symmetric
                        }
                }


/*
                multi3d<LatticeColorMatrix> f0;
                f0.resize(4,Nd,Nd);
                
                f0 = 0;
                
                for(int mu = 0; mu < Nd; mu++)
                {
                        for(int nu = mu+1; nu < Nd; nu++)
                        {
			for(int i = 0; i < 4; i++)
                        {
                        f0[i][nu][mu] = plane_plaq[i][nu][mu] - adj(plane_plaq[i][nu][mu]);
			f0[i][nu][mu] *= 1/8.0;
			f0[i][mu][nu] = -f0[i][nu][mu];
			}
                   	}
                }
		
		
		multi3d<LatticeColorMatrix> OP3;
		OP3.resize(5,4,4);
		OP3 = 0;
                for(int dir=0; dir<3; dir++)
		{
			for(int len=0; len<5; len++)
			{
				OP3 = fun_Operator3(len, dir, f0, u[dir]);
                        for(int t = 0; t < Layout::lattSize()[3]; t++)
                        {
				multi3d<Double> op_3;
				op_3.resize(5,4,4);
				op_3=0;
                                multi1d<int> tCoords;
                                tCoords.resize(Nd);
                                tCoords[3] = t;
                                for(int x = 0; x < Layout::lattSize()[0];x++)
                                  for(int y = 0; y < Layout::lattSize()[1];y++)
                                    for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                                        {
                                                tCoords[0] = x;
                                                tCoords[1] = y;
                                                tCoords[2] = Z;
						for(int k = 0; k<5;k++)
						for(int i = 0; i<4;i++)
						for(int j = 0; j<4;j++)
						{
                                                        op_3[k][i][j] += real(trace(peekSite(OP3[k][i][j], tCoords)));
						}


                                        }
				for(int k = 0; k<5;k++)
				for(int i = 0; i<4;i++)
				for(int j = 0; j<4;j++)
				{
                                	QDPIO::cout <<"Op0123   dir="<< dir << "  len="<< len <<"  t="<<  t <<"  subdir="<< k <<"  p1="<< i <<"  p2="<< j <<"  "<< op_3[k][i][j] <<std::endl;
				}
			}

			}
		}
*/		

/*

		LatticeColorMatrix Op0, Op1, Op2, Op3, Op4, Op5, Op6, Op7, Op8, Op9, Op10, Op11;
                
                
			Op0 = f0[0][dir][3]*f0[0][dir][3];
                        Op1 = f0[0][dir][3]*f0[1][dir][3];
                        Op2 = f0[0][dir][3]*f0[2][dir][3];
                        Op3 = f0[0][dir][3]*f0[3][dir][3];
			//Op4 = shift(f0[1][dir][3],FORWARD,dir)*shift(f0[1][dir][3],FORWARD,dir);
			//Op5 = f0[0][dir][3]*shift(f0[1][dir][3],FORWARD,dir);
			Op4 = f0[0][dir][3]*u[dir]*shift(u[dir],FORWARD,dir)*shift(shift(f0[0][dir][3],FORWARD,dir),FORWARD,dir)*adj(shift(u[dir],FORWARD,dir))*adj(u[dir]);
                        Op5 = f0[0][dir][3]*u[dir]*shift(u[dir],FORWARD,dir)*shift(shift(f0[1][dir][3],FORWARD,dir),FORWARD,dir)*adj(shift(u[dir],FORWARD,dir))*adj(u[dir]);
                        Op6 = f0[0][dir][3]*u[dir]*shift(u[dir],FORWARD,dir)*shift(shift(f0[2][dir][3],FORWARD,dir),FORWARD,dir)*adj(shift(u[dir],FORWARD,dir))*adj(u[dir]);
                        Op7 = f0[0][dir][3]*u[dir]*shift(u[dir],FORWARD,dir)*shift(shift(f0[3][dir][3],FORWARD,dir),FORWARD,dir)*adj(shift(u[dir],FORWARD,dir))*adj(u[dir]);
                        Op8 = f0[2][dir][3]*u[dir]*shift(u[dir],FORWARD,dir)*shift(shift(f0[0][dir][3],FORWARD,dir),FORWARD,dir)*adj(shift(u[dir],FORWARD,dir))*adj(u[dir]);
                        Op9 = f0[2][dir][3]*u[dir]*shift(u[dir],FORWARD,dir)*shift(shift(f0[1][dir][3],FORWARD,dir),FORWARD,dir)*adj(shift(u[dir],FORWARD,dir))*adj(u[dir]);
                        Op10 = f0[2][dir][3]*u[dir]*shift(u[dir],FORWARD,dir)*shift(shift(f0[2][dir][3],FORWARD,dir),FORWARD,dir)*adj(shift(u[dir],FORWARD,dir))*adj(u[dir]);
                        Op11 = f0[2][dir][3]*u[dir]*shift(u[dir],FORWARD,dir)*shift(shift(f0[3][dir][3],FORWARD,dir),FORWARD,dir)*adj(shift(u[dir],FORWARD,dir))*adj(u[dir]);





                        for(int t = 0; t < Layout::lattSize()[3]; t++)
                        {
                                Double op0, op1, op2, op3, op4, op5, op6, op7, op8, op9, op10, op11;
                                
                                op0=0;
                                op1=0;
                                op2=0;
                                op3=0;
                                op4=0;
                                op5=0;
                                op6=0;
                                op7=0;
				op8=0;
				op9=0;
				op10=0;
				op11=0;
                                
                                multi1d<int> tCoords;
                                tCoords.resize(Nd);
                                tCoords[3] = t;
                                for(int x = 0; x < Layout::lattSize()[0];x++)
                                  for(int y = 0; y < Layout::lattSize()[1];y++)
                                    for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                                        {
                                                tCoords[0] = x;
                                                tCoords[1] = y;
                                                tCoords[2] = Z;
                                                
                                                        op0 += real(trace(peekSite(Op0, tCoords)));
                                                        op1 += real(trace(peekSite(Op1, tCoords)));
                                                        op2 += real(trace(peekSite(Op2, tCoords)));
                                                        op3 += real(trace(peekSite(Op3, tCoords)));
                                                        op4 += real(trace(peekSite(Op4, tCoords)));
                                                        op5 += real(trace(peekSite(Op5, tCoords)));
                                                        op6 += real(trace(peekSite(Op6, tCoords)));
                                                        op7 += real(trace(peekSite(Op7, tCoords)));
                                                        op8 += real(trace(peekSite(Op8, tCoords)));
                                                        op9 += real(trace(peekSite(Op9, tCoords)));
                                                        op10 += real(trace(peekSite(Op10, tCoords)));
                                                        op11 += real(trace(peekSite(Op11, tCoords)));


                                                
                                        }
                                QDPIO::cout <<"Op0123   "<< dir << "    "<< t <<"  "<< op0 <<"  "<< op1 <<"  "<< op2 <<"  "<< op3 <<"  "<< op4<<"  "<< op5<<"  "<< op6<<"  "<< op7<<"  "<< op8<<"  "<< op9<<"  "<< op10<<"  "<< op11<<std::endl;

                        }//end of t loop

*/


		
/*
                for(int len = 0; len < mn; len++)
                {

                	ColorMatrix wl_forward, wl_backward, wl_tem;
                        multi1d<LatticeColorMatrix> Op;
                        Op.resize(6);
                        Op = 0;
                        for(int t = 0; t < Layout::lattSize()[3]; t++)
                        {
                                multi1d<Double> op;
                                op.resize(6);
                                op = 0;
                                multi1d<int> tCoords, tCoordsl;
				int t_tem;
                                tCoords.resize(Nd);
                                tCoordsl.resize(Nd);
                                tCoords[3] = t;
                                tCoordsl[3] = t;
                                for(int x = 0; x < Layout::lattSize()[0];x++)
                                  for(int y = 0; y < Layout::lattSize()[1];y++)
                                    for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                                        {
                                                tCoordsl[0] = x;
                                                tCoordsl[1] = y;
                                                tCoordsl[2] = Z;
                                                tCoords[0] = x;
                                                tCoords[1] = y;
                                                tCoords[2] = Z;
						wl_forward = 1;
						wl_backward = 1;
						
						for(int k = 0; k < len; k++)
						{
							wl_tem = wl_forward;
							wl_forward = wl_tem*peekSite(u[dir], tCoords);
							wl_tem = wl_backward;
							wl_backward = peekSite(adj(u[dir]), tCoords)*wl_tem;
                                                        t_tem = tCoords[dir];
                                                        tCoords[dir] = (t_tem+1)%(Layout::lattSize()[dir]);
						
						}
						
                                                        for(int i=0;i<4; i++)
                                                        {
								op[0] += real(trace(peekSite(F[3][i], tCoordsl)*wl_forward*peekSite(F[3][i], tCoords)*wl_backward));
                                                                op[1] += real(trace(peekSite(F[3][i], tCoordsl)*wl_forward*peekSite(F[dir][i], tCoords)*wl_backward));
								op[2] += real(trace(peekSite(F[dir][i], tCoordsl)*wl_forward*peekSite(F[dir][i], tCoords)*wl_backward));
                                                                op[3] += real(trace(peekSite(F[dir][i], tCoordsl)*wl_forward*peekSite(F[dir][i], tCoords)*wl_backward));
                                                                op[4] += real(trace(peekSite(F[3][i], tCoordsl)*wl_forward*peekSite(F[3][i], tCoords)*wl_backward));
								if(i==dir || i==3)
								{op[5] += real(trace(peekSite(F[dir][i], tCoordsl)*wl_forward*peekSite(F[dir][i], tCoords)*wl_backward));}
                                                        }
                                                        for(int i=0;i<4; i++)
                                                        {
								for(int j=0;j<i;j++)
								{
									op[0] -= 0.5*real(trace(peekSite(F[j][i], tCoordsl)*wl_forward*peekSite(F[j][i], tCoords)*wl_backward));
									op[2] -= 0.5*real(trace(peekSite(F[j][i], tCoordsl)*wl_forward*peekSite(F[j][i], tCoords)*wl_backward));
								}
							}

                                        }
                                QDPIO::cout <<"Op   "<< dir << "  " << len <<"  "<< t <<"  "<< op[0] <<"  "<< op[1] <<"  "<< op[2] <<"  "<< op[3] <<"  "<< op[4]<<"  "<< op[5]<< std::endl;

                        }//end of t loop
		}//end of len loop
                
*/

		
		// Obtain the F_21 at t=8, y=12, z=12 along x direction
                std::vector<Double> VecFmnre, VecFmnim;
                for(int x = 0; x < Layout::lattSize()[0]; x++)
                {
                        Double fmnre = 0, fmnim = 0;
                        multi1d<int> fCoords;
                        fCoords.resize(Nd);
                        fCoords[3] = 8;
                        fCoords[0] = x;
                        fCoords[1] = 12;
                        fCoords[2] = 12;
                        fmnre = real(trace(peekSite(F[2][1], fCoords)));
                        fmnim = imag(trace(peekSite(F[2][1], fCoords)));
                        VecFmnre.push_back(fmnre);
			VecFmnim.push_back(fmnim);
                }

		// Print the F_\mu\nu at t=8, y=12, z=12 along x direction
                for(int x = 0; x < Layout::lattSize()[0]; x++)
                {
                       // QDPIO::cout <<"F21   "<< x << "   " << VecFmnre.at(x) << "   " << VecFmnim.at(x) << std::endl;

                }

		// Obtain and print the sum(trace( of F_\mu\nu to check the anti-symmetric.
                for(int t = 0; t < Layout::lattSize()[3]; t++)
                {
                	multi2d<Double> F_sumre, F_sumim; //real part and imaginary part of the sum(trace(F_\mu\nu))
			F_sumre.resize(Nd,Nd);
			F_sumim.resize(Nd,Nd);
			F_sumre = 0;
                        F_sumim = 0;
                 	multi1d<int> tCoords;
                        tCoords.resize(Nd);
                        tCoords[3] = t;
                	for(int x = 0; x < Layout::lattSize()[0]; x++)
                	{
			 for(int y = 0; y < Layout::lattSize()[0]; y++)
                	 {
                          for(int Z = 0; Z < Layout::lattSize()[0]; Z++)
                          {
			  	tCoords[0] = x;
				tCoords[1] = y;
				tCoords[2] = Z;
				for(int mu = 0; mu < Nd; mu++)
                		{
                        	 for(int nu = 0; nu < Nd; nu++)
                        	 {
					F_sumre[mu][nu] += real(trace(peekSite(F[mu][nu], tCoords)));
					F_sumim[mu][nu] += imag(trace(peekSite(F[mu][nu], tCoords)));
				 }
				}	
			  }	
			 }
			}
                        for(int mu = 0; mu < Nd; mu++)
                        {
                        	for(int nu = 0; nu < Nd; nu++)
                        	{
					//QDPIO::cout <<"F   "<< t << "   "<< mu << "   " << nu << "   "  << F_sumre[mu][nu] << "   " << F_sumim[mu][nu] << std::endl;
				}
			}
          	}


		// loop over the wilson line length len
		for(int len = 0; len < mn; len++)
        	{
			multi1d<LatticeColorMatrix> Op;
			Op.resize(13);
			Op = 0;
			Op = fun_Operator(len, dir, F, u[dir]); //Obtain the operators




                        multi2d<LatticeColorMatrix> Opcomp;
                        Opcomp.resize(3,Nd);
                        Opcomp = 0;
                        Opcomp = fun_Opcomp(len, dir, F, u[dir]); //Obtain the operators components



			//Calculate and print out the real(trace of the operators
			for(int t = 0; t < Layout::lattSize()[3]; t++)
                        {
                                multi1d<Double> op, opi;
				op.resize(11);
				op = 0;
                                opi.resize(11);
                                opi = 0;


                                multi2d<Double> opcomp;
                                opcomp.resize(3,Nd);
                                opcomp = 0;


                                multi1d<int> tCoords;
                                tCoords.resize(Nd);
                                tCoords[3] = t;
                                for(int x = 0; x < Layout::lattSize()[0];x++)
                                  for(int y = 0; y < Layout::lattSize()[1];y++)
                                    for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                                        {
                                                tCoords[0] = x;
                                                tCoords[1] = y;
                                                tCoords[2] = Z;
						for(int i = 0; i < 11;i++)
						{
                                                	op[i] += real(trace(peekSite(Op[i], tCoords)));
							opi[i] += imag(trace(peekSite(Op[i], tCoords)));
						}

						for(int k = 0; k < 3; k++)
						{
                                                for(int l = 0; l < 4; l++)
                                                {
                                                        opcomp[k][l] += real(trace(peekSite(Opcomp[k][l], tCoords)));
                                                }
						}

					}
				/*
                		for(int l = 0; l < 4; l++)
                		{
                        		QDPIO::cout <<"Opcomp_zz   "<< dir << "  " << len <<"  "<< t <<"  "<< l <<"  "<< opcomp[0][l]<<std::endl;
                        		QDPIO::cout <<"Opcomp_tz   "<< dir << "  " << len <<"  "<< t <<"  "<< l <<"  "<< opcomp[1][l]<<std::endl;
                        		QDPIO::cout <<"Opcomp_tt   "<< dir << "  " << len <<"  "<< t <<"  "<< l <<"  "<< opcomp[2][l]<<std::endl;

                		}*/


//				QDPIO::cout <<"Op   "<< dir << "  " << len <<"  "<< t <<"  "<< op[0] <<"  "<< op[1] <<"  "<< op[2] <<"  "<< op[3] <<"  "<< op[4]<<"  "<< op[5]<< "  "<< op[6]<<"  "<< op[7]<<"  "<< op[8]<<"  "<< op[9]<<"  "<< op[10]<<"  "<< op[11]<<"  "<< op[12]<<std::endl;
				QDPIO::cout <<"O0ZF   "<< dir << "  " << len <<"  "<< t <<"  "<< op[0] <<"  "<< opi[0]<<std::endl;
                                QDPIO::cout <<"O1ZF   "<< dir << "  " << len <<"  "<< t <<"  "<< op[1] <<"  "<< opi[1]<<std::endl;
                                QDPIO::cout <<"O2ZF   "<< dir << "  " << len <<"  "<< t <<"  "<< op[2] <<"  "<< opi[2]<<std::endl;
                                QDPIO::cout <<"O3ZF   "<< dir << "  " << len <<"  "<< t <<"  "<< op[3] <<"  "<< opi[3]<<std::endl;
                                QDPIO::cout <<"OIB    "<< dir << "  " << len <<"  "<< t <<"  "<< op[4] <<"  "<< opi[4]<<std::endl;
                                QDPIO::cout <<"O1JZ   "<< dir << "  " << len <<"  "<< t <<"  "<< op[5] <<"  "<< opi[5]<<std::endl;
                                QDPIO::cout <<"O2JZ   "<< dir << "  " << len <<"  "<< t <<"  "<< op[6] <<"  "<< opi[6]<<std::endl;
                                QDPIO::cout <<"O3JZ   "<< dir << "  " << len <<"  "<< t <<"  "<< op[7] <<"  "<< opi[7]<<std::endl;
                                QDPIO::cout <<"O1hJZ  "<< dir << "  " << len <<"  "<< t <<"  "<< op[8] <<"  "<< opi[8]<<std::endl;
                                QDPIO::cout <<"O2hJZ  "<< dir << "  " << len <<"  "<< t <<"  "<< op[9] <<"  "<< opi[9]<<std::endl;
                                QDPIO::cout <<"O3hJZ  "<< dir << "  " << len <<"  "<< t <<"  "<< op[10] <<"  "<< opi[10]<<std::endl;
				
			}//end of t loop
		}//end of len loop






       /* 
                for(int len = 1; len < 2; len++)
                {
                        multi1d<LatticeColorMatrix> Op;
                        Op.resize(6);
                        Op = 0;
                        for(int t = 0; t < Layout::lattSize()[3]; t++)
                        {
                                multi1d<Double> op;
                                op.resize(6);
                                op = 0;
                                multi1d<int> tCoords, tCoords1, tCoords2;
                                tCoords.resize(Nd);
				tCoords1.resize(Nd);
				tCoords2.resize(Nd);
                                tCoords[3] = t;
				tCoords1[3] = t;
				tCoords2[3] = t;
                                for(int x = 0; x < Layout::lattSize()[0];x++)
                                  for(int y = 0; y < Layout::lattSize()[1];y++)
                                    for(int Z = 0; Z < Layout::lattSize()[2];Z++)
                                        {
                                                tCoords[0] = x;
                                                tCoords[1] = y;
                                                tCoords[2] = Z;
						tCoords1[0] = (x+1)%(Layout::lattSize()[0]);
                                                tCoords1[1] = y;
                                                tCoords1[2] = Z;
                                                tCoords2[0] = (x+2)%(Layout::lattSize()[0]);
                                                tCoords2[1] = y;
                                                tCoords2[2] = Z;
                                                
                					for(int i=0;i<4; i++)
                					{
                        					op[1] += real(trace(peekSite(F[3][i], tCoords)*peekSite(u[0], tCoords)*peekSite(F[0][i], tCoords1)*peekSite(adj(u[0]), tCoords)));
								op[2] += real(trace(peekSite(F[3][i], tCoords)*peekSite(u[0], tCoords)*peekSite(u[0], tCoords1)*peekSite(F[0][i], tCoords2)*peekSite(adj(u[0]), tCoords1)*peekSite(adj(u[0]), tCoords)));;
               				 		}
                                        }
                                QDPIO::cout <<"Op   "<< 0 << "  " << len <<"  "<< t <<"  "<< op[0] <<"  "<< op[1] <<"  "<< op[2] <<"  "<< op[3] <<"  "<< op[4]<<"  "<< op[5]<< std::endl;

                        }//end of t loop
                }//end of len loop
	*/



	}//end of dir loop


    } 

}
                   
