#include <cassert>  // for assert
#include <cmath>    // for std::sin
#include <iomanip>  // for std::ios::fixed, std::setiosflags, std::setprecision
#include <iostream> // for std::cout, std::endl
#include <utility>  // for std::make_pair, std::pair 

namespace{
    static auto constexpr MAXBUFSIZE = 32;

    template <typename T, bool check_positive_num = true>
    T input_parameter(std::string const & str){
        T param;
        while (true){
		    std::cout << str;
            std::cin >> param;

		    std::cin.clear();
		    std::cin.ignore(MAXBUFSIZE, '\n');

            // 初期条件として<=0の場合に再入力が要らない時は<type, false>
            if constexpr (check_positive_num){
                if (!std::cin.fail() && param > static_cast<T>(0)){ 
			        break;
		        }
            } 
            else{
                if (!std::cin.fail()){
                    break;
                }
            }
	    }

        return param;
    }

    std::pair<double, double> integ_interval();

    double Simpson(int n, double x0, double x1);

    double f(double x);
}

int main(){
    auto const [low_end, up_end] = integ_interval();
    auto const n = input_parameter<int> ("Enter the partition number\n");
    auto const res = Simpson(n, low_end, up_end);

    std::cout << std::setprecision(7) << std::setiosflags(std::ios::fixed);
    std::cout << "When the partition number is " << n << ", the numerical integral in ["
    << low_end << "," << up_end << "] of y = sin(x) is " << res << std::endl; 
}

namespace{
    std::pair<double, double> integ_interval(){
        double low_end, up_end;

        while (true){
            // 積分の下端と上端の入力
		    std::cout << "Enter endpoints of for the integration interval" << std::endl;
		    std::cin >> low_end >> up_end;
            
		    if (!std::cin.fail() && low_end < up_end){
			    break;
		    }

		    std::cin.clear();
		    std::cin.ignore(MAXBUFSIZE, '\n');
	    }
    
        return std::make_pair(low_end, up_end);
    }

    double Simpson(int n, double x0, double x1){
        // x0 >= x1の時, 強制終了
        assert(x0 < x1);

        // 分割幅
        auto const dh = (x1 - x0) / static_cast<double>(n); 
        // 面積の初期化
        double sum = 0.0; 

        // Simpson's rule
        for (auto i = 1; i <= 2 * n; i++) {
            auto const x = x0 + 0.5 * static_cast<double>(i) * dh;

            if (i == 0 || i == 2 * n){
                sum += f(x);
            }
            else if (i % 2 == 0){
                sum += 2 * f(x);
            }
            else{
                sum += 4 * f(x);
            }
        }

        return (dh / 6) * sum;
    }

    double f(double x){
        return std::sin(x);
    }
}
