#ifndef PAYOFF_H
#define PAYOFF_H


#include <algorithm> // This is needed for the std::max comparison function, used in the pay-off calculations

class PayOff {
  protected : 
    double _E; //Strike

  public :
    PayOff(); // Default (no parameter) constructor 
    virtual ~PayOff() {}; // Virtual destructor

    // // Overloaded () operator, turns the PayOff into an abstract function object
    // virtual double operator() (const double& S) const = 0;

    // virtual pure compute function  - fixed strike
    virtual double computeFixed(const double& mean) const = 0; 
    
    // virtual pure compute function - floated strike
    virtual double computeFloated(const double& mean, const double& S) const = 0;   
};

class PayOffCall : public PayOff { 

    public :
        PayOffCall ( const double& E) ; 
        virtual ~PayOffCall () {};
        
        // surcharge of the compute function - fixed strike
        virtual double computeFixed(const double& mean) const; 
          // surcharge of the compute function - floated strike
        virtual double computeFloated(const double& mean, const double& S) const; 
        
    

};


#endif // PAYOFF_H
