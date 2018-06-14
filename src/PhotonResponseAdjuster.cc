#include "PhotonResponseAdjuster.h"

PhotonResponseAdjuster aPhotonResponseAdjuster ;

PhotonResponseAdjuster::PhotonResponseAdjuster() : Processor("PhotonResponseAdjuster") {

  // modify processor description
  _description = "Adjusts Photons from filtered pandora pfos collection" ;

  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "Printing" ,
                              "Print certain messages"  ,
                              _printing,
                               (int)5 ) ;
  
  registerProcessorParameter( "EnergyScaleFactor",
			      "Energy Scale Factor",
			      _energyScaleFactor,
			      (double)1.0);
//matching parameters
  registerProcessorParameter( "SmearAngles",
				"option 0 or 1 to smear the angles using smearing model and MC photon direction",
				_smearAngles,
				(int) 0);

  registerProcessorParameter( "AngularSmearingModel",
				"Energy dependence choice: 0 - stochastic, 1 - energy independent",
				_angularSmearingModel,
				(int) 0);

  registerProcessorParameter( "dTheta",
			      "normally distributed theta std dev (in units of radians) for drawing random variates in new direction ",
			       _dTheta,
			       (double) 0.001);

  registerProcessorParameter( "dPhi",
				"normally distributed phi std dev (in units of radians) for drawing random variates in new direction ",
				_dPhi,
				(double) 0.001);

  registerProcessorParameter( "AllowedEnergyDeviation",
			      " allowed energy deviation of MC and REC photon in standard deviations",
			      _allowedEnergyDeviation,
			      (double)999.0);
 
  registerProcessorParameter( "AllowedAngularDeviation",
			      " allowed theta deviation of MC and REC photon in standard deviations",
			      _allowedAngularDeviation,
			      (double)3.14);

  std::string inputParticleCollectionName = "PandoraPFOs";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                             "InputParticleCollectionName" ,
                             "Input Particle Collection Name "  ,
                             _inputParticleCollectionName,
                              inputParticleCollectionName);

  registerInputCollection( LCIO::MCPARTICLE,
                           "MCParticleCollection" ,
                           "Name of the MCParticle input collection"  ,
                           _mcParticleCollectionName ,
                           std::string("MCParticle") ) ;

 // std::string outputParticleCollectionName = "FastReconstructedParticles";
  std::string outputParticleCollectionName = "AdjustedPhotons";
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                             "OutputParticleCollectionName" ,
			     "Output Particle Collection Name "  ,
                             _outputParticleCollectionName,
                             outputParticleCollectionName);
}

void PhotonResponseAdjuster::init() {

  streamlog_out(DEBUG) << "   init called  "
                       << std::endl ;

  // usually a good idea to
  printParameters() ;
  rng = new TRandom3();
  _nEvt = 0;
  nrejected =0;

}

void PhotonResponseAdjuster::processRunHeader( LCRunHeader* run) {
  streamlog_out(MESSAGE) << " processRunHeader "  << run->getRunNumber() << std::endl ;
  _nRun++ ;
}

bool PhotonResponseAdjuster::FindPFOs( LCEvent* evt ) {

  bool tf = false;
std::cout<<" here? "<<std::endl;
  // clear old vector
  _pfovec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){
 // std::cout<<"collection item: "<< *itname << " with " << _inputParticleCollectionName<<std::endl;    
    if(*itname==_inputParticleCollectionName){
//      std::cout<< "found matching collection name for photon" <<std::endl;
      LCCollection* col = evt->getCollection(*itname);
      unsigned int nelem = col->getNumberOfElements();
//      std::cout<<" number of elements "<<nelem<<std::endl;
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
        ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
        _pfovec.push_back(recoPart);
      }
    }
  }

  if(_printing>1)std::cout << "Find PFOs : " << tf << std::endl;

  return tf;
}
bool PhotonResponseAdjuster::FindMCParticles( LCEvent* evt ){
         bool tf = false;

  // clear old vector
	 _mcpartvec.clear();
	 _mcpartflags.clear();

         typedef const std::vector<std::string> StringVec ;
         StringVec* strVec = evt->getCollectionNames() ;
         for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){
               if(*itname==_mcParticleCollectionName){
                     LCCollection* col = evt->getCollection(*itname);
                     unsigned int nelem = col->getNumberOfElements();
                      tf = true;
                      for(unsigned int i=0;i<nelem;i++){
                          MCParticle* mcPart = dynamic_cast<MCParticle*>(col->getElementAt(i));
                          _mcpartvec.push_back(mcPart);
			  _mcpartflags.push_back(false);      // initialization
                      }
                }
          }
  
          if(_printing>1)std::cout << "Find MCParticles : " << tf << " size " << _mcpartvec.size()<<std::endl; 
  
          return tf;
 }
double* PhotonResponseAdjuster::resimulateDirection(TLorentzVector mcgamma){
    
//	double * PhotonSmearing::SmearAngles(const double* oldmom, double sigma_theta, double sigma_phi){
       // TVector3 oldVec(oldmom[0], oldmom[1], oldmom[2]);
       // double theta=oldVec.Theta();
       // do1uble phi=oldVec.Phi();
//current model (still being validated)
//dTheta is millrad deviation
//do something like normal distribution width is something like dTheta/sqrt(E_mc )

        double sigma1;
        double sigma2;

// Choose the angular smearing model
        if(_angularSmearingModel == 0){                     // Stochastic
           sigma1 = _dTheta/std::sqrt(mcgamma.E());
	   sigma2 = _dPhi/std::sqrt(mcgamma.E());
        }
        else{                                               // Energy Independent
           sigma1 = _dTheta;
           sigma2 = _dPhi;
        }

/*	double phi = mcgamma.Phi();
        double theta = mcgamma.Theta();
        double dx = _dTheta/std::sqrt(mcgamma.E());//temp stuff to see if this mode actually works
        double dy = _dPhi/std::sqrt(mcgamma.E());
        double r = 2000; // assume resolution in mm (based on scale of detector thingies) say r is approx 2 meters out
        //trying new spherical model based on resoultion from x,y 
        //eventually rename dtheta dphi to dx and dy
        double a1 = cos(phi)*cos(phi)*dx*dx - sin(phi)*sin(phi)*dy*dy;
        double b1 = ( pow(cos(phi),4) - pow(sin(phi),4))*r*r*cos(theta)*cos(theta) ;
        double sigma1 = std::sqrt(a1/b1); // this is for dtheta

        double a2 = sin(phi)*sin(phi)*dx*dx - cos(phi)*cos(phi)*dy*dy;
        double b2 = ( pow(sin(phi),4) - pow(cos(phi),4))*r*r*sin(theta)*sin(theta);
        double sigma2 = std::sqrt(a2/b2) ; //this is for dphi
*/
 //       double newTheta = rng->Gaus(mcgamma.Theta(), sigma1);
  //      double newPhi = rng->Gaus(mcgamma.Phi(), sigma2);
        double newTheta = rng->Gaus(mcgamma.Theta(), sigma1);
        double newPhi = rng->Gaus(mcgamma.Phi(), sigma2);

//	std::cout<<"e t p "<<mcgamma.E()<<" "<<mcgamma.Theta()<<" "<<mcgamma.Phi()<<std::endl;
//	std::cout<<"sigma 1 2 "<<sigma1<<" "<<sigma2<<std::endl;
	//new unit vector direction
        TVector3 newVec(sin(newTheta)*cos(newPhi),sin(newTheta)*sin(newPhi),cos(newTheta));
        double* mom = new double[3];
        mom[0]=newVec.X();
        mom[1]=newVec.Y();
        mom[2]=newVec.Z();

        return mom;

}

double PhotonResponseAdjuster::safeAcos(double x){
	if (x < -1.0) x = -1.0 ;
	else if (x > 1.0) x = 1.0 ;
	return acos (x) ;
 }
float* PhotonResponseAdjuster::getNewCovMatrix(int pdg, double Energy){
	
	float* cov = new float[6];
	float E = (float)Energy;
	//initialize the entire matrix to 0s
	for(int i=0; i<6; i++){
		cov[i]=0.0;
	}
	if(pdg == 22){	//covariant terms are 0
		cov[0]=(0.18*std::sqrt(E) );
      
	}else{
		cov[0]=(0.55*std::sqrt(E) );
	}
	cov[2]=(0.001/std::sqrt(E) );
	cov[5]=(0.001/std::sqrt(E) );
	return cov;
}

int PhotonResponseAdjuster::getCorrespondingMCParticleIndex(TLorentzVector rec){

// This needs to be reworked.    Graham.   In busy events it appears to currently pick up 
// the first available match, and then leaves no particle to match with for subsequent ReconstructedParticles.
// 
// Matching would probably be better if based simply on angular deviation.

        if(_mcpartvec.size() == 0) return -1;
        int closest_match_index=-1;
        double angular_residual=-1;
        double e_residual=0;
        double tempresidual1=-1;
        double tempresidual2=-1;
        TLorentzVector mc;

        for(unsigned int i=0; i<_mcpartvec.size(); i++){
                //compare angles
                if(_mcpartflags.at(i) == true) continue;

                mc.SetPxPyPzE(_mcpartvec[i]->getMomentum()[0],_mcpartvec[i]->getMomentum()[1],_mcpartvec[i]->getMomentum()[2],_mcpartvec[i]->getEnergy());

                tempresidual1 = (rec.Vect()).Angle(mc.Vect())/0.0015;             // Calculate angle in space between RP and MCP
                tempresidual2 = fabs( (rec.E() - mc.E())/(0.18*sqrt(mc.E())) );   // Calculate number of standard deviations

        	std::cout<<"residuals "<<tempresidual1<<" "<<tempresidual2<<std::endl;	 
                       if((closest_match_index==-1) &&
				(tempresidual2  <= _allowedEnergyDeviation) &&
				(tempresidual1 <= _allowedAngularDeviation)) {
                                closest_match_index = i;
                                angular_residual = tempresidual1;
                                e_residual = tempresidual2;
                        }
                   
                        double bestSoFar = angular_residual*angular_residual + e_residual*e_residual;
                        double currentOne = tempresidual1*tempresidual1 + tempresidual2*tempresidual2;

                        if(( currentOne < bestSoFar ) &&

				(tempresidual2 <= _allowedEnergyDeviation) &&
				(tempresidual1 <= _allowedAngularDeviation) ){

                                closest_match_index=i;
                                angular_residual = tempresidual1;
                                e_residual = tempresidual2;
                        }
        }
        if(closest_match_index != -1){
                mc.SetPxPyPzE(_mcpartvec[closest_match_index]->getMomentum()[0],_mcpartvec[closest_match_index]->getMomentum()[1],_mcpartvec[closest_match_index]->getMomentum()[2],_mcpartvec[closest_match_index]->getEnergy());
                if(_printing>3){
                        std::cout<<"MC Match: "<<std::endl;
                        std::cout<<"Reco (E,theta,phi) "<<rec.E()<<", "<<rec.Theta()<<", "<<rec.Phi()<<" "<<std::endl;
                        std::cout<<"MC   (E,theta,phi) "<<mc.E()<<", "<<mc.Theta()<<", "<<mc.Phi()<<" " <<std::endl;
			_mcpartflags.at(closest_match_index) = true;
                }
        }
        else{
                if(_printing>3){
                std::cout<<"Photon not matched "<<std::endl;
                std::cout<<"Reco (E,theta,phi)"<<rec.E()<<", "<<rec.Theta()<<", "<<rec.Phi()<<" "<<std::endl;
		std::cout<<"total # rejected photons "<< nrejected++ <<std::endl;
                }
        }
//particle is matched so flag it
//	_mcpartflags.at(closest_match_index) = true;	
//	std::cout<<"returning match index "<<closest_match_index<<std::endl;
        return closest_match_index;
}

void PhotonResponseAdjuster::processEvent( LCEvent * evt ) {
 std::cout<<"starting to process event"<<std::endl;
 std::cout << "======================================== event " << _nEvt << std::endl ;
  //FindMCParticles(evt);
   FindPFOs(evt);
   if(_smearAngles) FindMCParticles(evt);

  // Make a new vector of particles
  LCCollectionVec * calreccol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  //fastreccol->setSubset(true);
 calreccol->setSubset(true);

  streamlog_out(MESSAGE) << " start processing event " << std::endl;

	//set up new particle
	double oldE, newE;
	double* newmom = new double[3];
	const double* oldmom;
	ParticleIDImpl* newPDG = new ParticleIDImpl();
	//newPDG->setPDG(22);
	//newPDG->setLikelihood(1.0);



	for(unsigned int i=0; i<_pfovec.size(); i++){
		ReconstructedParticleImpl* calRecoPart = new ReconstructedParticleImpl();

		if(_pfovec.at(i)->getType() == 22){
			//set photon specific pdg stuff
			newPDG->setPDG(22);
			newPDG->setLikelihood(1.0);

			//calculate adjusted energy and momentum
			oldE = _pfovec.at(i)->getEnergy();
			newE = _energyScaleFactor*oldE; 
			oldmom = _pfovec[i]->getMomentum();

			//use old unit vector and multiply by new energy
			newmom[0] = oldmom[0]/oldE;
			newmom[1] = oldmom[1]/oldE;
			newmom[2] = oldmom[2]/oldE;
		
			newmom[0] = newmom[0]*newE;
			newmom[1] = newmom[1]*newE;
			newmom[2] = newmom[2]*newE;

			/************************
			fix the photon direction based on MC 				
			direction use the already semi calibrated 				
			photon for better matching

			This section was for addressing the errors `			
			photon angle reconstruction, the bug has
			since been fixed.  This function really should really have its own method
			*************************/
			if(_smearAngles){
				std::cout<<"Smearing photon direction using MC direction"<<std::endl;
				TLorentzVector gamma;
				gamma.SetPxPyPzE(newmom[0],newmom[1],newmom[2],newE);
				TLorentzVector mcgamma;
				int mcindex;
				double* newdirection;
				mcindex = getCorrespondingMCParticleIndex(gamma);
				if(mcindex != -1){
					mcgamma.SetPxPyPzE(_mcpartvec.at(mcindex)->getMomentum()[0], _mcpartvec.at(mcindex)->getMomentum()[1], _mcpartvec.at(mcindex)->getMomentum()[2], _mcpartvec.at(mcindex)->getEnergy() );	
					newdirection = resimulateDirection(mcgamma); 
//					std::cout<<"matched and recomputing"<<std::endl;
					newmom[0] = newdirection[0]*newE;
					newmom[1] = newdirection[1]*newE;
					newmom[2] = newdirection[2]*newE;		
					TLorentzVector v;
					v.SetPxPyPzE(newmom[0],newmom[1],newmom[2],newE);
					
				}
			}//end angle smearing	

			//set up the recoparticle and add it to the
			//collection
			calRecoPart->setMomentum(newmom);
			calRecoPart->setEnergy(newE);
			//give the reco part an E,theta,phi cov matrix
			calRecoPart->setCovMatrix(getNewCovMatrix(22, newE) );
			calRecoPart->setMass(0.0);
			calRecoPart->setCharge(0.0);
			calRecoPart->addParticleID(newPDG);
			calRecoPart->setParticleIDUsed(newPDG);
			calRecoPart->setType(22);

		        calreccol->addElement( calRecoPart );
			
		}//end if 22
		else{	//if its not a photon copy over old elements
			newPDG->setPDG(_pfovec.at(i)->getType());
			newPDG->setLikelihood(1.0);
			//manually copy over the momentum, the array
			// needs its own space on the heap
			oldmom = _pfovec.at(i)->getMomentum();
			for(int i=0; i<3; i++){
				newmom[i] = oldmom[i];
			}
			calRecoPart->setMomentum(newmom);
			calRecoPart->setEnergy(_pfovec.at(i)->getEnergy());
			calRecoPart->setMass(_pfovec.at(i)->getMass());
			calRecoPart->addParticleID(newPDG);
			calRecoPart->setParticleIDUsed(newPDG);
			calRecoPart->setType(_pfovec.at(i)->getType());
			
			//if its not a photon but a different neutral particle make a new covariance matrix
			 if(_pfovec.at(i)->getCharge() == 0){
			//if this is not a photon or track,
			//add a covariance matrix for E,theta,phi
				calRecoPart->setCovMatrix(getNewCovMatrix(_pfovec.at(i)->getType(), _pfovec.at(i)->getEnergy()) );
			}
			calreccol->addElement( calRecoPart );	
			

		}//end else other neutrals


			//print out the changes to see the changes
		if(_printing>1){
		
		 	TLorentzVector gold,gnew;
			gold.SetPxPyPzE(oldmom[0],oldmom[1],oldmom[2],oldE);
			gnew.SetPxPyPzE(newmom[0],newmom[1],newmom[2],newE);

			std::cout<<"Event No. :"<< _nEvt << std::endl;
			std::cout<<"Adjusted Particle Pdg: "<< calRecoPart->getType() << std::endl;
			std::cout<<"Old (Px,Py,Pz,E): "<<
			oldmom[0] <<" "<<
			oldmom[1] <<" "<<
			oldmom[2] <<" "<<
			oldE<<std::endl;

			std::cout<<"Old (E,Theta,Phi): "<< 				
			gold.E()<< " "<< 
			gold.Theta() << " " << 
			gold.Phi() << std::endl;

			std::cout<<"New (Px,Py,Pz,E): "<< 
			newmom[0] <<" "<<
			newmom[1]<<" "<<
			newmom[2]<<" "<<
			newE<<std::endl;

			std::cout<<"New (E,Theta,Phi): "<<
			gnew.E()<< " "<<
			gnew.Theta() << " "<<
			gnew.Phi() << std::endl; 
			}
			
		
		
	}//end pfovec loop	
 
  _nEvt++;

  // Add new collection to event
  evt->addCollection( calreccol , _outputParticleCollectionName.c_str() );


 
}

void PhotonResponseAdjuster::end(){

}






