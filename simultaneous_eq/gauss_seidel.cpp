// Gauss Seidel method

#include <cassert>  // for assert
#include <cmath>    // for std::fabs, std::pow, std::sqrt
#include <iostream> // for std::cerr, std::cin, std::cout
#include <optional> // for std::optional
#include <string>   // for std::string
#include <utility>  // for std::move
#include <vector>   // for std::vector

namespace{
    static auto constexpr MAXBUFSIZE = 32;
    static auto constexpr MAXLOOP = 100000;
    static auto const eps = std::pow(2.0, -50);

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

    double calc_norm(std::vector<double> const &afterx, std::vector<double> const &beforex);

    std::optional<std::vector<double>> Gauss_Seidel(std::vector<std::vector<double>> const &A, std::vector<double> const &b, double eps);

    double inner_product(std::vector<double> const &a, std::vector<double> const &b);
}

int main(){
    auto const size = input_parameter<int> ("Enter the size of a coefficient matrix\n");

    // Ax = b 
    std::vector<std::vector<double>> A(size, (std::vector<double>(size)));
    std::vector<double> b(size);
    
    std::cout << "Enter elements of the coefficient matrix" << std::endl;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            A[i][j] = input_parameter<double, false>("");
        }
    }

    std::cout << "Enter elements of a right side matrix" << std::endl;
    for (int i = 0; i < size; i++){
        b[i] = input_parameter<double, false>("");
    }

    // Gauss Seidel methodが成功して解が収束していれば, resに解が入る
    if (auto const res = Gauss_Seidel(A, b, eps); res){
        // resから中身を取り出す
        auto const x = *res;

        for (auto i = 0; i < size; i++){
            std::cout << x[i] << std::endl;
        }
    }
    else{  
        // resに中身が入っていなかったら、解が収束していないのでエラーメッセージ
        std::cerr << "Error! Gauss_Seidel method failed!" << std::endl;

        return -1; 
    }
}

namespace{
    double calc_norm(std::vector<double> const &afterx, std::vector<double> const &beforex){
        auto const size = afterx.size();
        
        // afterxとbeforexのサイズが違っていたら強制終了
        assert(size == beforex.size());

        auto norm = 0.0;

        for (auto i =0ULL; i < size; i++){
            norm += (afterx[i] - beforex[i]) * (afterx[i] - beforex[i]);
        }

        return std::sqrt(norm);
    }

    std::optional<std::vector<double>> Gauss_Seidel(std::vector<std::vector<double>> const &A, std::vector<double> const &b, double eps){
        auto const size = A.size();

        // vect Aとvect bのサイズが違っていたら強制終了
        assert(size == b.size());

        std::vector<double> x(size, 0.0);

        for (auto i = 0; i < MAXLOOP; i++){
            auto const beforex(x);

            for (auto j =0ULL; j < size; j++){
                x[j] = (b[j] - (inner_product(A[j], x) - A[j][j] * x[j])) / A[j][j];
            }

            // 残差vectorのnormがEPS未満になったら終了
            if (calc_norm(x, beforex) < eps){
                // 収束したので、解ベクトルxをstd::moveで右辺値参照にキャストし、std::optionalでラップして返す
                return std::make_optional(std::move(x));
            }
        }

        return std::nullopt;
    }

    double inner_product(std::vector<double> const &a, std::vector<double> const &b){
        auto value = 0.0;

        auto const size = a.size();
        
        // vect aとvect bのサイズが違っていたら強制終了
        assert(size == b.size());

        for (auto i = 0ULL; i < size; i++){
            value += a[i] * b[i];
        }

    return value;
    }
}
