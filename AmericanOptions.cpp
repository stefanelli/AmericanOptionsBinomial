#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

using namespace std;


class BinomialTree {
  /* Main class that defines the recombining binomial tree and the option valuation function */
private: 

  /* Model inputs */
  int n, option_type;
  double volatility, stock_value, strike, risk_free_rate, divident_yield, expiry;

  /* Calculation of model parameters */
  double step() {   //Time step
    return expiry / n;
  }

  double up_multiplier() {  // Multiplier for an upward movement of the stock price
    double time_step = step();
    return exp(volatility * sqrt(time_step));
  }

/* Utility functions */

vector<double> compare_max(const vector<double> vv1, const vector<double> vv2) {
// Performs a element-wise comaprison of two given vector using the "max" function
// returns a vector of element-wise max
  assert(vv1.size() == vv2.size());
  int n = vv1.size();
  vector<double> max_vector (n, 0); 
  for(int i = 0; i < n ; ++i){
    max_vector[i] = max(vv1[i], vv2[i]);
  }
  return max_vector;
}


  vector<double> stock_vector_build(int k){
// Builds a vector of stock values for a given depth "k" of the binomial tree.
// Returns a vector of stock values
    vector<double> stock (k +1, 0);
    double up = up_multiplier();
    for(int i = 0; i < k+1; ++i){
      stock[i] = get_stock_value() * pow(up, k-2*i);
      }
    return stock;
  }

vector<double> expected_disc_value(const vector<double> option_vector){
// Calculates the vector of expected discounted values for the level N-1 of the binomial tree
// Input: a vector of option prices of size N; Output a vector of option prices of size N-1
  int k = option_vector.size(); 
  vector<double> result_vector (k-1, -1);
  double time_step = step();
  double up = up_multiplier();
  double a = exp((risk_free_rate - divident_yield) * time_step);
  double down = 1.0 / up;
  double up_probability = (a - down)/(up - down);
  double discount_factor = exp(-risk_free_rate * time_step);
  for(int i = 0; i < k-1; ++i){
    result_vector[i] = discount_factor * (up_probability * option_vector[i] + (1 - up_probability) * option_vector[i+1]);
    }
  return result_vector;
}

vector<double> intrinsic_value(int k, double strike) {
// Calculates the intrinsic value of the option at any given level of the Binomial tree
// Returns a vector
    vector<double> stock_vector = stock_vector_build(k);
    vector<double> result_vector (k+1, strike);
    for(int i =0; i < k+1; ++i){
      result_vector[i] = (stock_vector[i] - strike) * option_type;
    }
  return result_vector;
}

vector<double> terminal_value(const int k, const double strike) {
// Calculates the terminal values for the option at expiry using the option payoff function
// returns a vector
    vector<double> stock_vector = stock_vector_build(k);
    vector<double> result_vector (k+1, strike);
    for(int i =0; i < k+1; ++i){
      result_vector[i] = max((stock_vector[i] - strike) * option_type, 0.0);
    }
  return result_vector;
}


public: 

// Get functions
  int get_depth() const {return n;}
  double get_stock_value() const {return stock_value;}
  double get_volatility() const {return volatility;}
  double get_rfr() const {return risk_free_rate;}
  double get_dividend() const {return divident_yield;}
  double get_strike() const {return strike;}
  
// Set functions
  void set_tree_depth(const int & m) {n = m;}
  void set_stock_value(const double & stock) {stock_value = stock;}
  void set_vol(const double & vol) { 
    assert(vol > 0);
    volatility = vol;
  }
  void set_type(const int & flag) {
    assert(flag == 1 || flag == -1);
    option_type = flag;
    }
  void set_rfr(const double & rfr) {risk_free_rate = rfr;}
  void set_dividend(const double & dividend) {divident_yield = dividend;}
  void set_expiry(const double & expiry_time) {expiry = expiry_time;}
  void set_strike(const double & level) {strike = level;}
  
/* Main valuation funtion */
double option_value(){ 
  /* 
  Applies the main pricing logic:
  1- gets the terminal value for the option payoff
  2- compares the Max between the intrinsic value of the option and the expected discounted value for all layers of the binomial tree
  3- returns the first node of the tree which is the price of the American option
  */
    double strike = get_strike();
    vector<double> option_vector = terminal_value(n, strike);
    vector<double> result_value = option_vector;
    for(int i = result_value.size() - 1; i > 0; i--){
      option_vector = expected_disc_value(option_vector);
      result_value = compare_max(intrinsic_value(i-1, strike), option_vector);
      //print_vector(option_vector);
    }
    assert(result_value.size() == 1);
    return result_value[0];
}

};

class Sensitivities {
  /*
  This class defines the sensitivities calculations using the "bump and reval" method and central difference.
  The greeks available are : stock Delta; vega; IR Delta and dividend Delta.
  */
  private:
    double bump = 1.0/10000;
  public:
    
    void set_bump(const double & size) {bump = size;}
    
    double central_difference(double result_up, double result_down) {
      return (result_up - result_down) / (2 * bump);
    }
    
    double delta_stock(BinomialTree tree) {
      double base_value = tree.get_stock_value();
      tree.set_stock_value(base_value + bump);
      double result_up = tree.option_value();
      tree.set_stock_value(base_value - bump);
      double result_down = tree.option_value();
      return central_difference(result_up, result_down);
    }
    
    double vega(BinomialTree tree) {
      double base_value = tree.get_volatility();
      tree.set_vol(base_value + bump);
      double result_up = tree.option_value();
      tree.set_vol(base_value - bump);
      double result_down = tree.option_value();
      return central_difference(result_up, result_down);
    }
    
    double div_delta(BinomialTree tree) {
      double base_value = tree.get_dividend();
      tree.set_dividend(base_value * (1 + bump));
      double result_up = tree.option_value();
      tree.set_dividend(base_value * (1 - bump));
      double result_down = tree.option_value();
      return central_difference(result_up, result_down);
    }
    
    double ir_delta(BinomialTree tree) {
      double base_value = tree.get_rfr();
      tree.set_rfr(base_value * (1 + bump));
      double result_up = tree.option_value();
      tree.set_rfr(base_value * (1 - bump));
      double result_down = tree.option_value();
      return central_difference(result_up, result_down);
    }
};

////////////////////////////////////////

int main() {
  
  cout << "Americal option pricer:\n" << endl;
  
  int n                  = 20 ;
  double volatility      = 0.1;
  double stock_value     = 900;
  double strike          = 1000;
  int option_type        = 1;  // 1 for Call, -1 for Put
  double risk_free_rate  = 0.05;
  double divident_yield  = 0.02;
  double expiry          = 2;
  
  cout << "Tree depth is: " << n << endl;
  cout << "Initial value for the stock is:  " << stock_value << endl;
  cout << "Strike price is:  " << strike << endl;
  if(option_type == 1)  {cout << "The option is a Call" << endl;}
  if(option_type == -1) {cout << "The option is a Put"  << endl;}
  cout << "The risk free rate is constant at: " << risk_free_rate*100 <<"%" << endl;
  cout << "The divident yield is: " << divident_yield*100 << "%" << endl;
  cout << "The time of expiry is " << expiry << " units of time" << endl;
  
  BinomialTree valuation_tree;
  
  valuation_tree.set_tree_depth(n);
  valuation_tree.set_vol(volatility);
  valuation_tree.set_dividend(divident_yield);
  valuation_tree.set_rfr(risk_free_rate);
  valuation_tree.set_type(option_type);
  valuation_tree.set_expiry(expiry);
  valuation_tree.set_stock_value(stock_value);
  valuation_tree.set_strike(strike);
  
  
  double result = valuation_tree.option_value();
  
  cout << "\nThe Option PV is: " << result << endl;
    
  Sensitivities greeks;
  
  double delta_1 = greeks.delta_stock(valuation_tree);
  cout << "\nThe stock Delta is: " << delta_1 << endl;
  
  double vega = greeks.vega(valuation_tree);
  cout << "\nThe Vega is: " << vega << endl;
  
  double ir_delta = greeks.ir_delta(valuation_tree);
  cout << "\nThe IR delta is: " << ir_delta << endl;
  
  double div_delta = greeks.div_delta(valuation_tree);
  cout << "\nThe Dividend delta is: " << div_delta << endl;

  return 0;

}
