#include "inline_Ob.h"
#include "chroma.h"
#include <math.h>
#include <complex>
#include <iostream>
#include <string>
#include "stdio.h"
#include <time.h>



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
	InlineObParams::InlineObParams() { frequency = 0; radius = 0; pmax = 0; srcs.resize(1); }
	
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
		//read(paramtop, "Multi_Src", srcs);

		if(paramtop.count("radius") == 1)
		    read(paramtop, "radius", radius);
		else
		    radius = 0;

		read(paramtop, "pmax", pmax);
		
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
/*
    LatticeReal Delta(int i, int j)
    {
	if(i==j)
		Delta=1.;
	else
		Delta=0.;
	return Delta;
    }

    Double lambda_tree(int rep, int alpha, int beta, int tau, multi1d<Double>)
    {
	if(rep==3)
		lambda_tree=(Nc*Nc-1)/sqrt(2)*(p[alpha]*p[alpha]+p[beta]*p[beta]);
	else
		lambda_tree=(Nc*Nc-1)/sqrt(2)*(2*(-cmplx(0,1))*Delta(beta,4)*p[alpha]*p[alpha]);
	return lambda_tree;	
    }
*/
    

    LatticeColorMatrix Ax(int Mu, LatticeColorMatrix u)
    {
        LatticeColorMatrix A_x;

        double g0=1.0;
        A_x=1/(2*g0)*((u-adj(u)-1/Nc*trace(u-adj(u))));

        return A_x;
    }




    LatticeReal pdotx(multi1d<int> p, multi1d<int> xsrc, multi1d<LatticeInteger> my_coord)
    {

	LatticeReal p_dot_x=0., p_shift=0.5;
        //const Real twopi = 6.283185307179586476925286;

        clock_t t1, t2, t3, t4;
/*
        t1=clock();

    	multi1d<LatticeInteger> my_coord(Nd);
    	for (int mu=0; mu < Nd; ++mu)
      	  my_coord[mu] = Layout::latticeCoordinate(mu);

	t2=clock();
	QDPIO::cout <<"t_coord   "<< (t2-t1) <<std::endl;
*/
        for(int mu=0;mu<Nd;mu++)
        {
		t1=clock();
                //if(mu==Mu)
                       // p_dot_x += LatticeReal(my_coord[mu]+p_shift)*twopi*Real(p[mu])/Layout::lattSize()[mu];
                //else
                p_dot_x += LatticeReal(my_coord[mu])*twopi*Real(p[mu])/Layout::lattSize()[mu];
		t2=clock();
		QDPIO::cout <<"t_pdotx   "<< (t2-t1) <<std::endl;
	
	}
        //t2=clock();
        //QDPIO::cout <<"t_pdotx   "<< t2-t1 <<std::endl;

	t3=clock();
        multi1d<int> tCoords;
        tCoords.resize(Nd);
        tCoords[0] = 1;
        tCoords[1] = 1;
        tCoords[2] = 1;
        tCoords[3] = 1;
        Real pdotx=0.;
        pdotx=peekSite(p_dot_x, tCoords);
        QDPIO::cout <<"p_dot_x   "<< pdotx <<std::endl;
	t4=clock();
	QDPIO::cout <<"pdotx print  "<< (t4-t3) <<std::endl;
	return p_dot_x;

    }

    Complex G3ptnew(int Mu, int Nu, LatticeColorMatrix u, LatticeColorMatrix u1, LatticeReal p_dot_x, LatticeColorMatrix Op, multi1d<int> p, int rmax, int rmax0)
    {
        Complex G_3pt=0;

        LatticeColorMatrix A_x, A_x0, A_x0c, umid, ux=u, ux1=u1;
        ColorMatrix  A_p, A_mp;
	LatticeColorMatrix Opc;

        const Real twopi = 6.283185307179586476925286;

        double g0=1.0;
        A_x=1/(2*g0)*((ux-adj(ux)-1/Nc*trace(ux-adj(ux))));
        A_x0=1/(2*g0)*((ux1-adj(ux1)-1/Nc*trace(ux1-adj(ux1))));

        multi1d<int> xcoords, xrcoords, xr0coords;
        xcoords.resize(Nd);
        xrcoords.resize(Nd);
        xr0coords.resize(Nd);

       	LatticeReal p_shift=0.5;
	LatticeReal p_dot_x0, p_dot_x1;
	LatticeInteger mask, mask0;
	LatticeComplex phase, phasem;

	
	p_dot_x0 = p_dot_x+LatticeReal(LatticeReal(p_shift)*twopi*Real(p[Mu])/Layout::lattSize()[Mu]);
        phase=cmplx(-sin(p_dot_x0),-cos(p_dot_x0));
	
        p_dot_x1 = p_dot_x+LatticeReal(LatticeReal(p_shift)*twopi*Real(p[Nu])/Layout::lattSize()[Nu]);
        phasem=cmplx(sin(p_dot_x1),-cos(p_dot_x1));

        for(int x = 0; x < Layout::lattSize()[0]; x++)
        {
                QDPIO::cout <<"x   "<< x <<std::endl;
                xcoords[0] = x;
                for(int y = 0; y < Layout::lattSize()[1]; y++)
                {
                        xcoords[1] = y;
                        for(int z = 0; z < Layout::lattSize()[2]; z++)
                        {
                                xcoords[2] = z;
                                for(int t = 0; t < Layout::lattSize()[3]; t++)
                                {
                                QDPIO::cout <<"t   "<< t <<std::endl;
                                xcoords[3] = t;
     				mask |= (Layout::latticeCoordinate(0)-x)*(Layout::latticeCoordinate(0)-x)+(Layout::latticeCoordinate(1)-y)*(Layout::latticeCoordinate(1)-y)+(Layout::latticeCoordinate(2)-z)*(Layout::latticeCoordinate(2)-z)+(Layout::latticeCoordinate(3)-t)*(Layout::latticeCoordinate(3)-t) <= rmax;
				Opc = where(mask, Op, LatticeColorMatrix(zero));
				mask0 |= (Layout::latticeCoordinate(0)-x)*(Layout::latticeCoordinate(0)-x)+(Layout::latticeCoordinate(1)-y)*(Layout::latticeCoordinate(1)-y)+(Layout::latticeCoordinate(2)-z)*(Layout::latticeCoordinate(2)-z)+(Layout::latticeCoordinate(3)-t)*(Layout::latticeCoordinate(3)-t) <= rmax0;
                                A_x0c = where(mask0, A_x0, LatticeColorMatrix(zero));
				G_3pt=G_3pt+trace(sum(Opc))*trace(peekSite(phase*A_x, xcoords)*sum(phasem*A_x0c));
                                }
                        }
                }
        }

    	return G_3pt;
    }



    Complex G3pt(int Mu, int Nu, LatticeColorMatrix u, LatticeColorMatrix u1, LatticeColorMatrix Op, multi1d<int> p, int rmax, int rmax0)
    {
        Complex G_3pt=0;

        LatticeColorMatrix A_x, A_x0, umid, ux=u, ux1=u1;
        ColorMatrix  A_p, A_mp;

        const Real twopi = 6.283185307179586476925286;

        double g0=1.0;
        A_x=1/(2*g0)*((ux-adj(ux)-1/Nc*trace(ux-adj(ux))));
        A_x0=1/(2*g0)*((ux1-adj(ux1)-1/Nc*trace(ux1-adj(ux1))));
        
        multi1d<int> xcoords, xrcoords, xr0coords;
        xcoords.resize(Nd);
	xrcoords.resize(Nd);
	xr0coords.resize(Nd);
	
	Real p_dot_x=0;
 
	for(int x = 0; x < Layout::lattSize()[0]; x++)
        {
		QDPIO::cout <<"x   "<< x <<std::endl;
         	xcoords[0] = x;
                for(int y = 0; y < Layout::lattSize()[1]; y++)
                {
                	xcoords[1] = y;
                	for(int z = 0; z < Layout::lattSize()[2]; z++)
                    	{
				xcoords[2] = z;
				for(int t = 0; t < Layout::lattSize()[3]; t++)
				{
				QDPIO::cout <<"t   "<< t <<std::endl;
				xcoords[3] = t;
				for(int dx = 0; dx*dx < 4*rmax*rmax; dx++)
				for(int dy = 0; dx*dx+dy*dy < 4*rmax*rmax; dy++)	
				for(int dz = 0; dx*dx+dy*dy+dz*dz < 4*rmax*rmax; dz++)
				for(int dt = 0; dx*dx+dy*dy+dz*dz+dt*dt < 4*rmax*rmax; dt++)
				{
				QDPIO::cout <<"dt   "<< dt <<std::endl;
                                for(int dx0 = 0; dx0*dx0 < 4*rmax0*rmax0; dx0++)
                                for(int dy0 = 0; dx0*dx0+dy0*dy0 < 4*rmax0*rmax0; dy0++)
                                for(int dz0 = 0; dx0*dx0+dy0*dy0+dz0*dz0 < 4*rmax0*rmax0; dz0++)
				for(int dt0 = 0; dx0*dx0+dy0*dy0+dz0*dz0+dt0*dt0 < 4*rmax0*rmax0; dt0++)
                                {                        
					//QDPIO::cout <<"dt0   "<< dt0 <<std::endl;			
					xrcoords[0] = x+dx-rmax;
					if(x+dx-rmax<0)
						xrcoords[0] = Layout::lattSize()[0]+x+dx-rmax;
					if(x+dx-rmax>Layout::lattSize()[0])
						xrcoords[0] = -Layout::lattSize()[0]+x+dx-rmax;
                                        xrcoords[1] = y+dy-rmax;
                                        if(y+dy-rmax<0)
                                                xrcoords[1] = Layout::lattSize()[1]+y+dy-rmax;
                                        if(y+dy-rmax>Layout::lattSize()[1])
                                                xrcoords[1] = -Layout::lattSize()[1]+y+dy-rmax;
                                        xrcoords[2] = z+dz-rmax;
                                        if(z+dz-rmax<0)
                                                xrcoords[2] = Layout::lattSize()[2]+z+dz-rmax;
                                        if(z+dz-rmax>Layout::lattSize()[2])
                                                xrcoords[2] = -Layout::lattSize()[2]+z+dz-rmax;
					xrcoords[3] = t+dt-rmax;
                                        if(t+dt-rmax<0)
                                                xrcoords[3] = Layout::lattSize()[3]+t+dt-rmax;
                                        if(t+dt-rmax>Layout::lattSize()[3])
                                                xrcoords[3] = -Layout::lattSize()[3]+t+dt-rmax;
                                        xr0coords[0] = x+dx0-rmax0;
                                        if(x+dx0-rmax0<0)
                                                xr0coords[0] = Layout::lattSize()[0]+x+dx0-rmax0;
                                        if(x+dx0-rmax0>Layout::lattSize()[0])
                                                xr0coords[0] = -Layout::lattSize()[0]+x+dx0-rmax0;
                                        xr0coords[1] = y+dy0-rmax0;
                                        if(y+dy0-rmax0<0)
                                                xr0coords[1] = Layout::lattSize()[1]+y+dy0-rmax0;
                                        if(y+dy0-rmax0>Layout::lattSize()[1])
                                                xr0coords[1] = -Layout::lattSize()[1]+y+dy0-rmax0;
                                        xr0coords[2] = z+dz0-rmax0;
                                        if(z+dz0-rmax0<0)
                                                xr0coords[2] = Layout::lattSize()[2]+z+dz0-rmax0;
                                        if(z+dz0-rmax0>Layout::lattSize()[2])
                                                xr0coords[2] = -Layout::lattSize()[2]+z+dz0-rmax0;
					xr0coords[3] = t+dt0-rmax0;
                                        if(t+dt0-rmax0<0)
                                                xr0coords[3] = Layout::lattSize()[3]+t+dt0-rmax0;
                                        if(t+dt0-rmax0>Layout::lattSize()[3])
                                                xr0coords[3] = -Layout::lattSize()[3]+t+dt0-rmax0;
					p_dot_x=twopi*(p[0]*dx0/Layout::lattSize()[0]+p[1]*dy0/Layout::lattSize()[1]+p[2]*dz0/Layout::lattSize()[2]+p[3]*dt0/Layout::lattSize()[3]);
					for(int mu=0; mu<Nd; mu++)
					{	
						if(mu==Mu)
							p_dot_x=p_dot_x-twopi*0.5*p[mu]/Layout::lattSize()[mu];
						if(mu==Nu)
                                                        p_dot_x=p_dot_x+twopi*0.5*p[mu]/Layout::lattSize()[mu];
					}

					G_3pt=G_3pt+cmplx(cos(p_dot_x),sin(p_dot_x))*peekSite(trace(Op), xrcoords)*trace(peekSite(A_x, xcoords)*peekSite(A_x0, xr0coords));
				}
				}
				}
			}
		}
	}
    return G_3pt;
    }


    Complex G2ptnew(int Mu, multi1d<int> p, LatticeColorMatrix A_x, int rmax)
    {
        Complex G_2pt;
        ColorMatrix  A_p, A_mp;

        clock_t t1, t2, t3, t4;

	Real p_dot_x=0;
	LatticeComplex AA=0;
	LatticeColorMatrix A_x1, Amid;

	A_x1=A_x;
        for(int dx = 0; dx*dx < rmax*rmax; dx++)
        {
                Amid=A_x1;
                A_x1=shift(Amid, BACKWARD, 0);
                A_x1=shift(Amid, BACKWARD, 1);
                A_x1=shift(Amid, BACKWARD, 2);
                //A_x1=shift(Amid, BACKWARD, 3);
        }

        //for(int dt = 0; dt*dt < rmax*rmax; dt++)
        for(int dt = 0; dt*dt < 33*33; dt++)
	{
                Amid=A_x1;
                A_x1=shift(Amid, BACKWARD, 3);
	}

        //for(int dx = -rmax+1; dx*dx < rmax*rmax; dx++)
        //for(int dx = -rmax+1; dx*dx < rmax*rmax; dx++)
	for(int dx = -rmax+1; dx < rmax-1; dx++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 0);
        //for(int dy = -rmax+1; dx*dx+dy*dy < rmax*rmax; dy++)
        //for(int dy = -rmax+1; dy*dy < rmax*rmax; dy++)
        for(int dy = -rmax+1; dy < rmax-1; dy++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 1);
        //for(int dz = -rmax+1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz++)
        //for(int dz = -rmax+1; dz*dz < rmax*rmax; dz++)
        for(int dz = -rmax+1; dz < rmax-1; dz++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 2);
        //for(int dt = -rmax+1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt++)
        //for(int dt = -rmax+1; dt < rmax; dt++)
        for(int dt = -33+1; dt < 32; dt++)
	//for(int dt = -rmax+1; dt*dt < rmax*rmax; dt++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 3);
        p_dot_x=twopi*(1.0*p[0]*dx/Layout::lattSize()[0]+1.0*p[1]*dy/Layout::lattSize()[1]+1.0*p[2]*dz/Layout::lattSize()[2]+1.0*p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x*A_x1);
        }
        }
        }
        }

	G_2pt=sum(AA);
	return G_2pt;
/*
	A_x1=A_x;
	p_dot_x=0;
	AA=cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
	
	for(int dx = 1; dx*dx < rmax*rmax; dx++)
	{
		Amid=A_x1;
		A_x1=shift(Amid, FORWARD, 0);
		p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0];
		AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
	}
	
	A_x1=A_x;
        for(int dx = -1; dx*dx < rmax*rmax; dx--)
        {
                Amid=A_x1;
                A_x1=shift(Amid, BACKWARD, 0);
                p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0];
                AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }

	A_x1=A_x;
        for(int dy = 1; dy*dy < rmax*rmax; dy++)
        {
                Amid=A_x1;
                A_x1=shift(Amid, FORWARD, 1);
                p_dot_x=twopi*(p[1]*dy/Layout::lattSize()[1];
                AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }

        A_x1=A_x;
        for(int dy = -1; dy*dy < rmax*rmax; dy--)
        {
                Amid=A_x1;
                A_x1=shift(Amid, BACKWARD, 1);
                p_dot_x=twopi*(p[1]*dy/Layout::lattSize()[1];
                AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }

        A_x1=A_x;
        for(int dz = 1; dz*dz < rmax*rmax; dz++)
        {
                Amid=A_x1;
                A_x1=shift(Amid, FORWARD, 2);
                p_dot_x=twopi*(p[2]*dz/Layout::lattSize()[2];
                AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }

        A_x1=A_x;
        for(int dz = -1; dz*dz < rmax*rmax; dz--)
        {
                Amid=A_x1;
                A_x1=shift(Amid, BACKWARD, 2);
                p_dot_x=twopi*(p[2]*dz/Layout::lattSize()[2];
                AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }

        A_x1=A_x;
        for(int dt = 1; dt*dt < rmax*rmax; dt++)
        {
                Amid=A_x1;
                A_x1=shift(Amid, FORWARD, 3);
                p_dot_x=twopi*(p[3]*dz/Layout::lattSize()[3];
                AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }

        A_x1=A_x;
        for(int dt = -1; dt*dt < rmax*rmax; dt--)
        {
                Amid=A_x1;
                A_x1=shift(Amid, BACKWARD, 3);
                p_dot_x=twopi*(p[3]*dz/Layout::lattSize()[3];
                AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }


	A_x1=A_x;
        for(int dx = 1; dx*dx < rmax*rmax; dx++)
	{
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 0);
        for(int dy = 1; dx*dx+dy*dy < rmax*rmax; dy++)
	{
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 1);
        for(int dz = 1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz++)
	{
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 2);
        for(int dt = 1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt++)
        {
	Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 3);
	p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
	AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
	}
	}
	}
	}

        A_x1=A_x;
        for(int dx = 1; dx*dx < rmax*rmax; dx++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 0);
        for(int dy = 1; dx*dx+dy*dy < rmax*rmax; dy++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 1);
        for(int dz = 1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 2);
        for(int dt = -1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = 1; dx*dx < rmax*rmax; dx++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 0);
        for(int dy = 1; dx*dx+dy*dy < rmax*rmax; dy++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 1);
        for(int dz = -1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 2);
        for(int dt = 1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = 1; dx*dx < rmax*rmax; dx++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 0);
        for(int dy = 1; dx*dx+dy*dy < rmax*rmax; dy++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 1);
        for(int dz = -1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 2);
        for(int dt = -1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = 1; dx*dx < rmax*rmax; dx++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 0);
        for(int dy = -1; dx*dx+dy*dy < rmax*rmax; dy--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 1);
        for(int dz = 1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 2);
        for(int dt = 1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = 1; dx*dx < rmax*rmax; dx++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 0);
        for(int dy = -1; dx*dx+dy*dy < rmax*rmax; dy--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 1);
        for(int dz = 1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 2);
        for(int dt = -1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = 1; dx*dx < rmax*rmax; dx++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 0);
        for(int dy = -1; dx*dx+dy*dy < rmax*rmax; dy--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 1);
        for(int dz = -1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 2);
        for(int dt = 1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = 1; dx*dx < rmax*rmax; dx++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 0);
        for(int dy = -1; dx*dx+dy*dy < rmax*rmax; dy--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 1);
        for(int dz = -1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 2);
        for(int dt = -1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }


        A_x1=A_x;
        for(int dx = -1; dx*dx < rmax*rmax; dx--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 0);
        for(int dy = 1; dx*dx+dy*dy < rmax*rmax; dy++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 1);
        for(int dz = 1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 2);
        for(int dt = 1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = -1; dx*dx < rmax*rmax; dx--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 0);
        for(int dy = 1; dx*dx+dy*dy < rmax*rmax; dy++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 1);
        for(int dz = 1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 2);
        for(int dt = -1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = -1; dx*dx < rmax*rmax; dx--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 0);
        for(int dy = 1; dx*dx+dy*dy < rmax*rmax; dy++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 1);
        for(int dz = -1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 2);
        for(int dt = 1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = -1; dx*dx < rmax*rmax; dx--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 0);
        for(int dy = 1; dx*dx+dy*dy < rmax*rmax; dy++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 1);
        for(int dz = -1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 2);
        for(int dt = -1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = -1; dx*dx < rmax*rmax; dx--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 0);
        for(int dy = -1; dx*dx+dy*dy < rmax*rmax; dy--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 1);
        for(int dz = 1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 2);
        for(int dt = 1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = -1; dx*dx < rmax*rmax; dx--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 0);
        for(int dy = -1; dx*dx+dy*dy < rmax*rmax; dy--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 1);
        for(int dz = 1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 2);
        for(int dt = -1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = -1; dx*dx < rmax*rmax; dx--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 0);
        for(int dy = -1; dx*dx+dy*dy < rmax*rmax; dy--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 1);
        for(int dz = -1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 2);
        for(int dt = 1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt++)
        {
        Amid=A_x1;
        A_x1=shift(Amid, FORWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }

        A_x1=A_x;
        for(int dx = -1; dx*dx < rmax*rmax; dx--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 0);
        for(int dy = -1; dx*dx+dy*dy < rmax*rmax; dy--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 1);
        for(int dz = -1; dx*dx+dy*dy+dz*dz < rmax*rmax; dz--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 2);
        for(int dt = -1; dx*dx+dy*dy+dz*dz+dt*dt < rmax*rmax; dt--)
        {
        Amid=A_x1;
        A_x1=shift(Amid, BACKWARD, 3);
        p_dot_x=twopi*(p[0]*dx/Layout::lattSize()[0]+p[1]*dy/Layout::lattSize()[1]+p[2]*dz/Layout::lattSize()[2]+p[3]*dt/Layout::lattSize()[3]);
        AA=AA+cmplx(cos(p_dot_x),sin(p_dot_x))*trace(A_x1*A_x);
        }
        }
        }
        }
*/	
    }


    ColorMatrix  Ap(int Mu, LatticeReal p_dot_x, multi1d<int> p, LatticeColorMatrix A_x)
    {
        Complex G_2pt;
        ColorMatrix  A_p, A_mp;
         
	clock_t t1, t2, t3, t4;
	
	LatticeComplex phase;
        t1=clock();	
	LatticeReal p_shift=0.5;	
	p_dot_x += LatticeReal(p_shift)*twopi*Real(p[Mu])/Layout::lattSize()[Mu];
        phase=cmplx(-sin(p_dot_x),-cos(p_dot_x));
        A_p=sum(phase*A_x);

        t3=clock();
        QDPIO::cout <<"t_Ap   "<< t3-t1 <<std::endl;


	LatticeColorMatrix AA, A_x1, Amid;

/*        A_x1=A_x;
        t1=clock();

	t3=clock();
	for(int i=1; i<100; i++)
	{
        Amid=A_x1;
	}
	t4=clock();
        QDPIO::cout <<"Amid   "<< t4-t3 <<std::endl;

	t3=clock();
	for(int i=1; i<100; i++)
        {
        A_x1=shift(Amid, FORWARD, Mu);
        }
	t4=clock();
	QDPIO::cout <<"t_shift   "<< t4-t3 <<std::endl;

	t3=clock();
	for(int i=1; i<100; i++)
        {
	AA=A_x1*sum(A_x)*A_x;
        }
        t4=clock();
        QDPIO::cout <<"t_AA1   "<< t4-t3 <<std::endl;

        t3=clock();
	for(int i=1; i<100; i++)
        {
        AA=A_x1*(sum(A_x)*A_x);
        }
        t4=clock();
        QDPIO::cout <<"t_AA1   "<< t4-t3 <<std::endl;


        t3=clock();
	for(int i=1; i<100; i++)
        {
        AA=sum(A_x)*A_x;
        }
        t4=clock();
        QDPIO::cout <<"t_AA3   "<< t4-t3 <<std::endl;

        t3=clock();
	for(int i=1; i<100; i++)
        {
        AA=LatticeColorMatrix(sum(A_x))*A_x;
        }
        t4=clock();
        QDPIO::cout <<"t_AA3   "<< t4-t3 <<std::endl;

        t3=clock();
        for(int i=1; i<100; i++)
        {
	AA=A_x1*A_x;
        }
        t4=clock();
        QDPIO::cout <<"t_AA2   "<< t4-t3 <<std::endl;
	QDPIO::cout <<"t_AA2   "<< (t4-t3)*100000 <<std::endl;
*/
        t2=clock();
        QDPIO::cout <<"t_sumx   "<< t2-t1 <<std::endl;

	return A_p;
    }

    Complex G2pt(int Mu, int Nu, LatticeReal p_dot_x, LatticeReal p_dot_x1, LatticeColorMatrix u, LatticeColorMatrix u1,  multi1d<int> p, multi1d<int> xsrc)
    {
	Complex G_2pt;
	//xsrc[mu]?????
        LatticeColorMatrix AA, A_x, A_x1, Amid, umid, ux=u, ux1=u1;
	ColorMatrix  A_p, A_mp;
        //A_x=1/(2*cmplx(0,1)*g0)*((u[tau]-adj(u[tau])-1/Nc*trace(u[tau]-adj(u))));
 
        //LatticeReal p_dot_x=0.,  p_dot_x1=0., p_shift=0.5;

        //const Real twopi = 6.283185307179586476925286;
	LatticeComplex phase, phasem;
	
	clock_t t1, t2, t3, t4;

	t1=clock();
/*
        for(int mu=0;mu<Nd;mu++)
        {
		if(mu==Mu)
			p_dot_x += LatticeReal(Layout::latticeCoordinate(mu)+p_shift-xsrc[mu])*twopi*p[mu]/Layout::lattSize()[mu];
		else
			p_dot_x += LatticeReal(Layout::latticeCoordinate(mu)-xsrc[mu])*twopi*p[mu]/Layout::lattSize()[mu];
                if(mu==Nu)
                        p_dot_x1 += LatticeReal(Layout::latticeCoordinate(mu)+p_shift-xsrc[mu])*twopi*p[mu]/Layout::lattSize()[mu];
		else
        		p_dot_x1 += LatticeReal(Layout::latticeCoordinate(mu)-xsrc[mu])*twopi*p[mu]/Layout::lattSize()[mu];
		//umid=ux;
		//ux=field(mu, xsrc[mu], umid);
                //umid=ux1;
                //ux1=field(mu, xsrc[mu], umid);

        }
        //phase=cmplx(cos(p_dot_x),-sin(p_dot_x));
	//phasem=cmplx(cos(p_dot_x),sin(p_dot_x));

	t2=clock();
	QDPIO::cout <<"t_pdotx   "<< t2-t1 <<std::endl;

	multi1d<int> tCoords;
        tCoords.resize(Nd);
	tCoords[0] = 1;
	tCoords[1] = 1;
	tCoords[2] = 1;
	tCoords[3] = 1;
	Real pdotx=0.;
	pdotx=peekSite(p_dot_x, tCoords);	
	QDPIO::cout <<"p_dot_x   "<< pdotx <<std::endl;	
*/
	double g0=1.0;
	phase=cmplx(-sin(p_dot_x),-cos(p_dot_x));
	phasem=cmplx(sin(p_dot_x1),-cos(p_dot_x1));		
	A_x=1/(2*g0)*((ux-adj(ux)-1/Nc*trace(ux-adj(ux))));
        A_p=sum(phase*A_x);
	A_x=1/(2*g0)*((ux1-adj(ux1)-1/Nc*trace(ux1-adj(ux1))));
        A_mp=sum(phasem*A_x);

	t3=clock();
        QDPIO::cout <<"t_FT   "<< t3-t1 <<std::endl;

	A_x1=A_x;
	t1=clock();
	Amid=A_x1;
	A_x1=shift(Amid, FORWARD, Mu);
	AA=A_x1*sum(A_x)*A_x;
	t2=clock();
	QDPIO::cout <<"t_sumx   "<< t2-t1 <<std::endl;


	QDPIO::cout <<"A_p   "<< Mu << "  "<< Nu << "  "<< trace(A_p) <<std::endl;
	QDPIO::cout <<"A_mp   "<< Mu << "  "<< Nu << "  "<< trace(A_mp) <<std::endl;
	
                        for(int i = 0; i < Nc; i++)
                        for(int j = 0; j < Nc; j++)
                        {
				QDPIO::cout <<"A_p   "<< Mu << "  "<< Nu << "  "<< i << "  "<< j << "  "<< A_p.elem().elem().elem(i,j).real() << "  "<< A_p.elem().elem().elem(i,j).imag() <<std::endl;
				QDPIO::cout <<"A_mp   "<< Mu << "  "<< Nu << "  "<< i << "  "<< j << "  "<< A_mp.elem().elem().elem(i,j).real() << "  "<< A_mp.elem().elem().elem(i,j).imag() <<std::endl;
                        }

	
        G_2pt=trace(A_p*A_mp);
        return G_2pt;
    }

/*

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

*/
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

        multi1d<int> tCoords;
        tCoords.resize(Nd);
        tCoords[0] = 12;
        tCoords[1] = 12;
        tCoords[2] = 12;
        tCoords[3] = 32;
        Complex FF=0.;

	if(len==0)
	{
        for(int mu = 0; mu < Nd; mu++)
	for(int nu = 0; nu < Nd; nu++)
        {
                for(int Mu = 0; Mu < Nd; Mu++)
		for(int Nu = 0; Nu < Nd; Nu++)
                {
		
        		//FF=trace(peekSite(F[mu][nu]*F[Mu][Nu], tCoords));
			FF=trace(sum(F[mu][nu]*F[Mu][Nu]));
			QDPIO::cout <<"FF   "<< len << "  "<< mu << "  "<< nu << "  "<<  Mu << "  "<< Nu << "  "<< real(FF) << "  "<< imag(FF) <<std::endl;			
                }

        }
	}
	else
	{
        for(int mu = 0; mu < Nd; mu++)
        for(int nu = 0; nu < Nd; nu++)
        {
                for(int Mu = 0; Mu < Nd; Mu++)
                for(int Nu = 0; Nu < Nd; Nu++)
                {

                        FF=trace(sum(F[mu][nu]*un*Fn[Mu][Nu]*adj(un)));
                        QDPIO::cout <<"FF   "<< len << "  "<< mu << "  "<< nu << "  "<<  Mu << "  "<< Nu << "  "<< real(FF) << "  "<< imag(FF) <<std::endl;
                }

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
                        F[mu][nu]= -F[nu][mu]; //anti-symmetric
                        }
                }


	//Ax(int Mu, LatticeColorMatrix u)
	LatticeColorMatrix A_x;
	for(int mu = 3; mu < Nd; mu++)
        {
	    for(int len = 0; len < 11; len++)	
	    {	
		multi1d<LatticeColorMatrix> Op;
                Op.resize(13);
                Op = 0;
                Op = fun_Operator(len, mu, F, u[mu]); //Obtain the operators


        	QDPIO::cout <<"OZF0   "<< mu  << "  "<< trace(sum(Op[0])) <<std::endl;
		QDPIO::cout <<"O2JZ   "<< mu  << "  "<< trace(sum(Op[6])) <<std::endl;
	    }
/*
                        for(int i = 0; i < Nc; i++)
                        for(int j = 0; j < Nc; j++)
                        {
				QDPIO::cout <<"OZF0   "<< mu << "  "<< i << "  "<< j << "  "<< sum(Op[0]).elem().elem().elem(i,j).real() << "  "<< sum(Op[0]).elem().elem().elem(i,j).imag() <<std::endl;                                                                
				QDPIO::cout <<"O2JZ   "<< mu << "  "<< i << "  "<< j << "  "<< sum(Op[6]).elem().elem().elem(i,j).real() << "  "<< sum(Op[6]).elem().elem().elem(i,j).imag() <<std::endl;
			}


		A_x=Ax(mu,u[mu]);
		multi1d<int> tCoords;
        	tCoords.resize(Nd);
		for(int x = 11; x < 14; x++)
		for(int y = 11; y < 14; y++)
		for(int z = 11; z < 14; z++)
		for(int t = 31; t < 34; t++)
		{
        		tCoords[0] = x;
        		tCoords[1] = y;
        		tCoords[2] = z;
        		tCoords[3] = t;
        		Real ax=0., axi=0.;
        		ax=real(trace(peekSite(A_x, tCoords)));
			axi=imag(trace(peekSite(A_x, tCoords)));
			for(int i = 0; i < Nc; i++)
			for(int j = 0; j < Nc; j++)
			{
				QDPIO::cout <<"Ax   "<< mu << "  "<< x << "  " << y <<"  "<< z <<"  "<< t <<"  "<< i <<"  "<< j <<"  "<< peekSite(A_x, tCoords).elem().elem().elem(i,j) <<std::endl;
			}
        		QDPIO::cout <<"A_x   "<< mu << "  "<< x << "  " << y <<"  "<< z <<"  "<< t <<"  "<< ax <<"  "<< axi <<std::endl;
		}
*/	
	}

	//Double G_2pt(u, multi2d<Double> p, multi2d<Double> xsrc)
	
	multi1d<int> p; 
	multi1d<int> xsrc;
	p.resize(Nd);
	xsrc.resize(Nd);
	p=0;
	xsrc=0;
	Complex GL2pt, GL3pt;
	
	clock_t t1, t2;
	const Real twopi = 6.283185307179586476925286;
	
        multi1d<LatticeInteger> my_coord(Nd);
        for (int mu=0; mu < Nd; ++mu)
          my_coord[mu] = Layout::latticeCoordinate(mu);

	int pmax=params.pmax;
	QDPIO::cout <<"pmax   "<< pmax <<std::endl;
	multi4d<LatticeReal> p_dot_x;
	p_dot_x.resize(2*pmax+1,2*pmax+1,2*pmax+1,2*pmax+1);
	//for(int mu = 0; mu < Nd; mu++)
	for(int i = -pmax; i < pmax+1; i++)
	for(int j = -pmax; j < pmax+1; j++)
	for(int k = -pmax; k < pmax+1; k++)
	for(int l = -pmax; l < pmax+1; l++)
	{
                p[0]=i;
                p[1]=j;
                p[2]=k;
                p[3]=l;
                xsrc[0]=0;
                xsrc[1]=0;
                xsrc[2]=0;
                xsrc[3]=0;
        	t1=clock();
		p_dot_x[i+pmax][j+pmax][k+pmax][l+pmax]=pdotx(p, xsrc, my_coord);
		t2=clock();
		QDPIO::cout <<"p_dot_x_p   "<< (t2-t1) <<std::endl;
		QDPIO::cout <<"p_dot_x_p   "<< (double)(t2-t1)/ CLOCKS_PER_SEC <<std::endl;
	}
	
	multi1d<LatticeColorMatrix> Op;	
	Op = fun_Operator(0, 0, F, u[0]);


        double g0=1.0;
	//LatticeColorMatrix A_x;
        multi5d<ColorMatrix> A_p;
        A_p.resize(Nd,2*pmax+1,2*pmax+1,2*pmax+1,2*pmax+1);

	for(int mu = 0; mu < Nd; mu++)
	{
        A_x=1.0/(2*g0)*((u[mu]-adj(u[mu])-1/Nc*trace(u[mu]-adj(u[mu]))));
        for(int i = -pmax; i < pmax+1; i++)
                for(int j = -pmax; j < pmax+1; j++)
                        for(int k = -pmax; k < pmax+1; k++)
                        for(int l = -pmax; l < pmax+1; l++)
                        {
                                p[0]=i;
                                p[1]=j;
                                p[2]=k;
                                p[3]=l;
                                xsrc[0]=0;
                                xsrc[1]=0;
                                xsrc[2]=0;
                                xsrc[3]=0;
                                t1=clock();
				A_p[mu][i+pmax][j+pmax][k+pmax][l+pmax]=Ap(mu, p_dot_x[i+pmax][j+pmax][k+pmax][l+pmax], p, A_x);
	                        t2=clock();
                                QDPIO::cout <<"time_A_p   "<< (t2-t1) <<std::endl;
                                QDPIO::cout <<"time_A_p   "<< (double)(t2-t1)/ CLOCKS_PER_SEC <<std::endl;
			
			}
	}

	for(int i = -pmax; i < pmax+1; i++)
        	for(int j = -pmax; j < pmax+1; j++)
			for(int k = -pmax; k < pmax+1; k++)
			for(int l = -pmax; l < pmax+1; l++)
			{
				p[0]=1;
				p[1]=0;
				p[2]=-2;
				p[3]=3;
				xsrc[0]=0;
				xsrc[1]=0;
				xsrc[2]=0;
				xsrc[3]=0;
				for(int mu = 0; mu < Nd; mu++)
				//for(int nu = 0; nu < Nd; nu++)
				{
					int nu;
					nu=mu;
					t1=clock();
					GL2pt=trace(A_p[mu][i+pmax][j+pmax][k+pmax][l+pmax]*A_p[nu][-i+pmax][-j+pmax][-k+pmax][-l+pmax]);
					//GL2pt=G2pt(mu, nu, p_dot_x[mu][i+pmax][j+pmax][k+pmax][l+pmax], p_dot_x[nu][i+pmax][j+pmax][k+pmax][l+pmax], u[mu], u[nu], p, xsrc);
					QDPIO::cout <<"G2pt   "<< mu << "  " << nu << "  "<< i << "  " << j <<"  "<< k <<"  "<< l <<"  "<< real(GL2pt) << "  " << imag(GL2pt) <<std::endl;
                                       	t2=clock();
                			QDPIO::cout <<"time_G2pt   "<< (t2-t1) <<std::endl;
                			QDPIO::cout <<"time_G2pt   "<< (double)(t2-t1)/ CLOCKS_PER_SEC <<std::endl;


					A_x=1.0/(2*g0)*((u[mu]-adj(u[mu])-1/Nc*trace(u[mu]-adj(u[mu]))));

					for(int r = 13; r < 14; r++)
					{
                                        t1=clock();
                                        GL2pt=G2ptnew(mu, p, A_x, r);
                                        QDPIO::cout <<"G2ptr   "  << r << "  " << mu << "  " << nu << "  "<< i << "  " << j <<"  "<< k <<"  "<< l <<"  "<< real(GL2pt) << "  " << imag(GL2pt) <<std::endl;
                                        t2=clock();
                                        QDPIO::cout <<"time_G2ptr   " << r << "  "<< (double)(t2-t1)/ CLOCKS_PER_SEC <<std::endl;
					}
/*
                                        t1=clock();
                                        GL2pt=G2ptnew(mu, p, A_x, 1);
                                        QDPIO::cout <<"G2ptr1   "<< mu << "  " << nu << "  "<< i << "  " << j <<"  "<< k <<"  "<< l <<"  "<< real(GL2pt) << "  " << imag(GL2pt) <<std::endl;
                                        t2=clock();
                                        QDPIO::cout <<"time_G2ptr1   "<< (double)(t2-t1)/ CLOCKS_PER_SEC <<std::endl;

                                        t1=clock();
                                        GL2pt=G2ptnew(mu, p, A_x, 2);
                                        QDPIO::cout <<"G2ptr2   "<< mu << "  " << nu << "  "<< i << "  " << j <<"  "<< k <<"  "<< l <<"  "<< real(GL2pt) << "  " << imag(GL2pt) <<std::endl;
                                        t2=clock();
                                        QDPIO::cout <<"time_G2ptr2   "<< (double)(t2-t1)/ CLOCKS_PER_SEC <<std::endl;

                                        t1=clock();
                                        GL2pt=G2ptnew(mu, p, A_x, 3);
                                        QDPIO::cout <<"G2ptr3   "<< mu << "  " << nu << "  "<< i << "  " << j <<"  "<< k <<"  "<< l <<"  "<< real(GL2pt) << "  " << imag(GL2pt) <<std::endl;
                                        t2=clock();
                                        QDPIO::cout <<"time_G2ptr3   "<< (double)(t2-t1)/ CLOCKS_PER_SEC <<std::endl;
*/

					//GL2pt=G2pt(mu, nu, p_dot_x[mu][i+pmax][j+pmax][k+pmax][l+pmax], p_dot_x[nu][i+pmax][j+pmax][k+pmax][l+pmax], u[mu], u[nu], p, xsrc);
					//QDPIO::cout <<"G2pt   "<< mu << "  " << nu << "  "<< i << "  " << j <<"  "<< k <<"  "<< l <<"  "<< real(GL2pt) << "  " << imag(GL2pt) <<std::endl;
					//t1=clock();
					//GL3pt=G3ptnew(mu, nu, u[mu], u[nu], p_dot_x[i+pmax][j+pmax][k+pmax][l+pmax], Op[0], p, 64, 64);
                                        //t2=clock();
                                        //QDPIO::cout <<"time_G3pt   "<< (t2-t1) <<std::endl;
                                        //QDPIO::cout <<"time_G3pt   "<< (double)(t2-t1)/ CLOCKS_PER_SEC <<std::endl;
					//QDPIO::cout <<"G3pt   "<< mu << "  " << nu << "  "<< i << "  " << j <<"  "<< k <<"  "<< l <<"  "<< real(GL3pt) << "  " << imag(GL3pt) <<std::endl;

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


	}//end of dir loop


    } 

}
                   
