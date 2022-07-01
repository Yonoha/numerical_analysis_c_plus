// Gauss Elimination

#include <cmath>    // for std::fabs, std::pow
#include <iostream> // for std::cerr, std::cin, std::cout
#include <optional> // for std::optional
#include <string>   // for std::string
#include <vector>   // for std::vector
#include <utility>  // for std::move
#include <boost/format.hpp>          // for boost::format
#include <boost/multi_array.hpp>     // for boost::multi_array
#include <boost/range/algorithm.hpp> // for boost::range::copy

namespace{
    using dmatrix = boost::multi_array<double, 2>;
    using dvector = std::vector<double>;

    static auto constexpr MAXBUFSIZE = 32;
    static auto const eps = pow(2, -50);

    template <typename T, bool check_positive_num = true>
    T input_parameter(std::string const &str){
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

    void display_result(dvector const &x);

    // 後退代入
    void Gauss_backward(dmatrix &A);

    // Gaussの消去法
    std::optional<dvector> Gauss_elimination(dmatrix &A, double eps);

    // 掃き出し法
    bool Gauss_forward(dmatrix &A, double eps);

    // pivot選択
    void Gauss_pivot(dmatrix &A, int k);

    dvector input_row(int size);

    dvector make_result(dmatrix const &A);    
}

int main(){
    auto const size = input_parameter<int> ("Enter the size of a expansion coefficient matrix n-rows and (n+1)-colowns\n");

    dmatrix A(boost::extents[size][size + 1]);
    
    std::cout << "Enter elements of the expansion coefficient matrix" << std::endl;
    for (int i = 0; i < size; i++){
        auto const row = input_row(size + 1);
        boost::range::copy(row, A[i].begin());
    }

    if (auto const res = Gauss_elimination(A, eps); res){
        auto const x(*res);
        display_result(x);
    }
    else{
        std::cerr << "Error! Gaussian Elimination failed" << std::endl;

        return -1;
    }
}

namespace{
    void display_result(dvector const &x){
        auto const size = static_cast<int>(x.size());
        for (int i = 0; i < size; i++){
            std::cout << "x[" << i << boost::format("] = %.7f\n") % x[i];
        }
    }

    void Gauss_backward(dmatrix &A){
        auto const size = static_cast<int>(A.size());

        A[size - 1][size] /= A[size - 1][size - 1];
        for (int i = size - 2; i >= 0; i--){
            auto ax = 0.0;
            for (int j = i + 1; j < size; j++){
                ax += A[i][j] * A[j][size];
            }

            A[i][size] = (A[i][size] - ax) / A[i][i];
        }
    }

    std::optional<dvector> Gauss_elimination(dmatrix &A, double eps){
        if (!Gauss_forward(A, eps)){
            return std::nullopt;
        }
        
        Gauss_backward(A);

        return std::make_optional(make_result(A));
    }

    bool Gauss_forward(dmatrix &A, double eps){
        auto const size = static_cast<int>(A.size());

        for (int k = 0; k < size - 1; k++) {
            Gauss_pivot(A, k);

            for (int i = k + 1; i < size; i++){
                for (int j = k; j < size + 1; j++){
                    A[i][j] -= A[i][k] * (A[k][j] / A[k][k]);
                    auto tmp = A[i][j];
                    auto a = 0;
                }
            }
        }

        if (std::fabs(A[size - 1][size - 1]) < eps){
            return false;
        }

        return true;
    }

    void Gauss_pivot(dmatrix &A, int k){
        auto const size = static_cast<int>(A.size());

        auto pivot_ren = k;
        auto pivot = std::fabs(A[k][k]);
        for (int i = k + 1; i < size; i++){
            if (std::fabs(A[i][k]) > pivot){
                pivot_ren = i;
                pivot = std::fabs(A[i][k]);
            }
        }

        if (k != pivot_ren){
            dvector tmp(A[k].begin(), A[k].end());
            boost::range::copy(A[pivot_ren], A[k].begin());
            boost::range::copy(tmp, A[pivot_ren].begin());
        }
    }

    dvector input_row(int size){
        dvector row(size);

        auto success = true;
        do{
            auto i = 0;
            for (; i < size; i++) {
                std::cin >> row[i];

                if (std::cin.fail()) {
                    success = false;
                    break;
                }
            }
            if (i == size) {
                success = true;
            }

            std::cin.clear();
            std::cin.ignore(MAXBUFSIZE, '\n');
        } while (!success);

        return row;
    }

    dvector make_result(dmatrix const &A){
        auto const size = static_cast<int>(A.size());

        dvector x(size);

        for (int i = 0; i < size; i++){
            x[i] = A[i][size];
        }

        return x;
    }
}
