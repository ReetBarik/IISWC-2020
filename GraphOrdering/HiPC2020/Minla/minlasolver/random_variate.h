#include <LEDA/core/array.h>

extern double gen_prob();

class random_variate {
private:
  array<int> I;
  array<double> D;
  bool is_double;
public:
  int generate()
  {
    if(is_double == true)
      {
        int i;
        double sum = 0;
        for(i=0; i<D.size(); i++)
          sum+=D[i];
        double p = sum * gen_prob();
        //        cerr << "sum=" << sum << ", " << p << endl;
        double prev = 0;
        sum = 0;
  
        for(i=0; i<D.size(); i++)
          {
            sum += D[i];
            //            cerr << prev << " --- " << p << " --- " << sum << endl;
            if((p>=prev)&&(p<=sum))
              {
                return i;
              }
            prev = sum;
          }
        return i-1;
      }
    else
      {
        
      }
  }
  
  random_variate(array<int> R)
  {
    I = R;
    is_double = false;
  }
  
  random_variate(array<double> R)
  {
    D = R;
    is_double = true;
  }
};
