// 1dim freefall by 4th Runge Kutta method

#include <iomanip>  // for std::ios::fixed, std::setiosflags, std::setprecision
#include <iostream> // for std::cout
#include <string>   // for std::string

namespace{
    static auto constexpr MAXBUFSIZE = 32;
    static auto constexpr c = 0.1;
    static auto constexpr g = 9.80665;

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

    // motion equation
    double f(double v);

    void Runge_Kutta(int n, double tmax, double v0);
}

int main(){
    // 分割数, 最大時間, 初期条件の入力
    auto const n = input_parameter<int> ("Enter the partition number n\n");
    auto const tmax = input_parameter<double> ("Enter the maximum time t_max[s]\n");
    auto const v0 = input_parameter<double, false> ("Enter the initial velocity v_0[m/s]\n");
    
    Runge_Kutta(n, tmax, v0);
}

namespace{
    double f(double v){
        return g - c * v;
    }

    void Runge_Kutta(int n, double tmax, double v0){
        // 分割幅
        auto const dt = tmax / static_cast<double>(n);
        // vの初期化
        auto v = v0;

        std::cout << std::setprecision(7) << std::setiosflags(std::ios::fixed);
        for (auto i = 0; i <= n; i++){
            auto const t = static_cast<double>(i) * dt;

            std::cout << "when the time is " << t << "[s], velocity is " << v << "[m/s]" << std::endl;

            // Runge Kutta method
            auto const k1 = f(v);
            auto const k2 = f(v + k1 * dt * 0.5);
            auto const k3 = f(v + k2 * dt * 0.5);
            auto const k4 = f(v + k3 * dt);
            v += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);            
        }
    }
}
