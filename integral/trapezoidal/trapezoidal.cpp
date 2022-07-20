#include <cassert>  // for assert
#include <iomanip>  // for std::ios::fixed, std::setiosflags, std::setprecision
#include <iostream> // for std::cout, std::endl
#include <utility>  // for std::make_pair, std::pair 

namespace{
    static auto constexpr MAXBUFSIZE = 32;

    std::pair<double, double> integ_interval();

    int partition_num();

    double trapezoidal_integral(int n, double x0, double x1);

    double y(double x);
}

int main(){
    auto const [low_end, up_end] = integ_interval();
    auto const n = partition_num();
    auto const res = trapezoidal_integral(n, low_end, up_end);

    std::cout << std::setprecision(7) << std::setiosflags(std::ios::fixed);
    std::cout << "When the partition number is " << n << ", the numerical integral in ["
    << low_end << "," << up_end << "] of y = x^2 is " << res << std::endl; 
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

    int partition_num(){
        int n;

        while (true){
            // 分割数の入力
		    std::cout << "Enter the partition number" << std::endl;
            std::cin >> n;

		    if (!std::cin.fail()) {
			    break;
		    }
            
		    std::cin.clear();
		    std::cin.ignore(MAXBUFSIZE, '\n');
	    }

        return n;
    }

    double trapezoidal_integral(int n, double x0, double x1){
        // x0 >= x1の時, 強制終了
        assert(x0 < x1);

        // 分割幅
        auto const dh = (x1 - x0) / static_cast<double>(n); 
        // 面積の初期化
        double sum = 0.0; 

        // trapezoidal integral
        for (auto i = 1; i < n; i++) {
            auto const x = x0 + static_cast<double>(i) * dh;
            sum += y(x);
        }

        return dh * ((y(x0) + y(x1)) / 2.0 + sum);
    }

    double y(double x){
        return x * x;
    }
}
