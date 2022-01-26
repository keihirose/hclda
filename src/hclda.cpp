// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace arma;  // use the Armadillo library for matrix computations
using namespace Rcpp;



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List creasteSS_Rcpp(
                    IntegerVector lab,
                    mat datx)
{
    
    int n = size(datx)[0];  //データのサンプル数
    int p = size(datx)[1];  //ベクトルの次元数
    CharacterVector lab_levels = lab.attr("levels");  //ラベル
    int G = lab_levels.length();  //群の数
    
    // n_g[g,1]でg群の要素数を表す
    uvec n_g(G, fill::zeros);
    for(int i = 0; i < n; ++i){
        n_g(lab(i)-1) += 1;
    }
    
    mat bar_x = mat(p, G, fill::zeros);
    for(int i = 0; i < G; ++i){
        for(int j = 0; j < n; ++j){
            if(lab(j) == i + 1){
                bar_x.col(i) = bar_x.col(i) + datx.row(j).t() / (double) n_g(lab(j)-1);
            }
        }
    }
    //colnames(bar_x) = levels(train_lab)
    vec ones_n(n, fill::ones);
    vec barx = datx.t() * ones_n / (double) n;
    
    //# 分散共分散行列と群間・郡内変動行列
    std::vector<mat> S_x(G);
    std::vector<mat> xjxjT(G);
    for(int i = 0; i < G; ++i){
        S_x[i].zeros(p, p);
        xjxjT[i].zeros(p, p);
    }
    for(int i = 0; i < G; ++i){
        for(int j = 0; j < n; ++j){
            if(lab(j) == i + 1){
                S_x[i] = S_x[i] + (datx.row(j).t() - bar_x.col(i)) * (datx.row(j).t() - bar_x.col(i)).t() / ((double) n_g(lab(j)-1)-1);
                xjxjT[i] = xjxjT[i] + (datx.row(j).t()) * datx.row(j);//  ここ要チェック！！！ / ((double) n_g(lab(i)-1)-1.0);
            }
        }
    }
    
    mat B = zeros(p,p);
    mat W = zeros(p,p);
    for(int i = 0; i < G; ++i){
        B = B + n_g(i) * (bar_x.col(i) - barx) * (bar_x.col(i) - barx).t();
        W = W + (n_g(i) - 1) * S_x[i];
    }
    //#以下，要チェック！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    //B <- B/n
    //W <- W/n
    mat tmp_pG = mat(p, G, fill::zeros);
    for(int i = 0; i < p; ++i){
        for(int j = 0; j < G; ++j){
            tmp_pG(i,j) = n_g(j);
        }
    }
    tmp_pG = tmp_pG % bar_x;
    mat sum_nj_xjbar_xjbar_T = tmp_pG * bar_x.t();
    
    
    Rcpp::Environment base("package:base");
    Rcpp::Function vector = base["vector"];
    Rcpp::Function is_element = base["is.element"];
    Rcpp::List S_x_list = vector(Rcpp::Named("length")=G, Rcpp::Named("mode")="list");
    Rcpp::List xjxjT_list = vector(Rcpp::Named("length")=G, Rcpp::Named("mode")="list");
    for(int i = 0; i < G; ++i){
        S_x_list[i] = S_x[i];
        xjxjT_list[i] = xjxjT[i];
    }
    
    
    return Rcpp::List::create(Rcpp::Named("G")=G,
                              Rcpp::Named("n.g")=n_g,
                              Rcpp::Named("bar.x")=bar_x,
                              Rcpp::Named("S.x")=S_x_list,
                              Rcpp::Named("xjxjT")=xjxjT_list,
                              Rcpp::Named("barx")=barx,
                              Rcpp::Named("B")=B,
                              Rcpp::Named("W")=W,
                              Rcpp::Named("sum_nj_xjbar_xjbar_T")=sum_nj_xjbar_xjbar_T);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat computeY_Rcpp(
                  int N,
                  int G,
                  int p,
                  double r,
                  mat W,
                  mat B,
                  mat bar_x,
                  vec barx,
                  uvec n_g)
{
    
    //# 固有値・固有ベクトル
    int D = std::min(p, G-1);
    vec eigval;
    mat eigvec;
    mat diagp = eye(p,p);
    eig_sym(eigval, eigvec, 0.5*(W+W.t()) + r * diagp);  // 固有値は小さい順であることに注意
    //eig_sym(eigval, eigvec, W + r * diagp);  // 固有値は小さい順であることに注意、上記は数値誤差を解消するため
    mat Winv_half = eigvec * diagmat(pow(abs(eigval), -0.5)) * eigvec.t();
    eig_sym(eigval, eigvec, 0.5* Winv_half * (B+B.t()) * Winv_half.t());
    //eig_sym(eigval, eigvec, Winv_half * B * Winv_half);
    // //eig_sym(eigval, eigvec, Winv_half * B_loo * Winv_half);
    mat hat_T = Winv_half * eigvec; //射影ベクトル A[,d]でd番目の固有ベクトル
    vec Val = eigval;
    
    //### 提案手法 ハット行列 ###
    //# Yを定義する
    mat Y = mat(N,D,fill::zeros);
    mat eps = mat(G,D,fill::zeros);
    int start = 0;
    int end = 0;
    for(int j = 0; j < G; ++j){
        for(int d = 0; d < D; ++d){
            //eps(j,d) = eps(j,d) + ( 1/((double) N * fabs(Val(p-d-1))) ) * n_g(k) * dot((bar_x.col(j) - bar_x.col(k)), hat_T.col(p-d-1));
            eps(j,d) = ( 1/(fabs(Val(p-d-1)))) * dot(bar_x.col(j) - barx,  hat_T.col(p-d-1));
        }
    }
    
    for(int j = 0; j < G; ++j){
        end = end + n_g(j);
        for(int k = start; k < end; k++){
            Y.row(k) = eps.row(j);
        }
        start = start + n_g(j);
    }
    
    //return Rcpp::List::create(Val, eigvec, Y, eps);
    return Y;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List creasteSS_fast_Rcpp(
                         IntegerVector lab,
                         mat datx,
                         double r)
{
    
    int n = size(datx)[0];  //データのサンプル数
    int p = size(datx)[1];  //ベクトルの次元数
    CharacterVector lab_levels = lab.attr("levels");  //ラベル
    int G = lab_levels.length();  //群の数
    
    // n_g[g,1]でg群の要素数を表す
    uvec n_g(G, fill::zeros);
    for(int i = 0; i < n; ++i){
        n_g(lab(i)-1) += 1;
    }
    
    mat bar_x = mat(p, G, fill::zeros);
    for(int i = 0; i < G; ++i){
        for(int j = 0; j < n; ++j){
            if(lab(j) == i + 1){
                bar_x.col(i) = bar_x.col(i) + datx.row(j).t() / (double) n_g(lab(j)-1);
            }
        }
    }
    //colnames(bar_x) = levels(train_lab)
    vec ones_n(n, fill::ones);
    vec barx = datx.t() * ones_n / (double) n;
    
    //# 分散共分散行列と群間・郡内変動行列
    std::vector<mat> S_x(G);
    std::vector<mat> xjxjT(G);
    for(int i = 0; i < G; ++i){
        S_x[i].zeros(p, p);
    }
    for(int i = 0; i < G; ++i){
        for(int j = 0; j < n; ++j){
            if(lab(j) == i + 1){
                S_x[i] = S_x[i] + (datx.row(j).t() - bar_x.col(i)) * (datx.row(j).t() - bar_x.col(i)).t() / ((double) n_g(lab(j)-1)-1);
            }
        }
    }
    
    mat B = zeros(p,p);
    mat W = zeros(p,p);
    for(int i = 0; i < G; ++i){
        B = B + n_g(i) * (bar_x.col(i) - barx) * (bar_x.col(i) - barx).t();
        W = W + (n_g(i) - 1) * S_x[i];
    }
    //#以下，要チェック！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    B = 1/(double)n * B;
    W = 1/(double)n * W;
    
    //最初の方のを追加, p>>nを想定する
    mat X_til(n, p+1, fill::ones);
    X_til.cols(1,p) = datx;
    mat Del = eye(p+1, p+1);
    Del = r * Del;
    Del(0,0)=0.0;
    mat XTX = X_til.t() * X_til;
    mat XTX_Del = XTX + Del;
    mat C_inv = inv(XTX_Del);
    mat ci = C_inv * X_til.t();
    mat Y = computeY_Rcpp(n, G, p, r, W, B, bar_x, barx, n_g);
    
    mat sum_y_tilx = X_til.t() * Y;
    mat alpha = C_inv * sum_y_tilx;
    mat hat_Y = X_til * alpha;
    vec sum_hatY2 = pow(hat_Y,2).t() * ones_n;
    mat sum_hY_H = X_til * (C_inv * (XTX * alpha));
    vec sum_H2(n, fill::zeros);
    vec diagH(n, fill::zeros);
    int D = std::min(p, G-1); //廣瀬変更！
    mat sum_hatY_j(G,D,fill::zeros); //廣瀬変更！
    mat sum_H_j(G,n,fill::zeros); //廣瀬変更！
    rowvec tmp;
    for(int j = 0; j < n; ++j){
        tmp = X_til.row(j) * ci;
        diagH(j) = tmp(j);
        for(int i = 0; i < G; ++i){
            if(lab(j) == i + 1){
                sum_hatY_j.row(i) = sum_hatY_j.row(i) + hat_Y.row(j);
                sum_H_j.row(i) = sum_H_j.row(i) + tmp;
            }
        }
        tmp = pow(tmp,2);
        sum_H2(j) = sum(tmp);
    }
    
    return Rcpp::List::create(Rcpp::Named("diagH")=diagH,
                              Rcpp::Named("hat_Y")=hat_Y,
                              Rcpp::Named("Y")=Y,
                              Rcpp::Named("alpha")=alpha,
                              Rcpp::Named("ci")=ci,
                              Rcpp::Named("sum_hatY2")=sum_hatY2,
                              Rcpp::Named("sum_hY_H")=sum_hY_H,
                              Rcpp::Named("sum_H2")=sum_H2,
                              Rcpp::Named("sum_hatY_j")=sum_hatY_j,
                              Rcpp::Named("sum_H_j")=sum_H_j,
                              Rcpp::Named("bar_x")=bar_x);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List loo_CLD_Rcpp(mat test_dat,
                  IntegerVector test_lab,
                  mat train_dat,
                  IntegerVector train_lab,
                  double r,
                  int N_press,
                  mat bar_x,
                  vec barx,
                  List xjxjT,
                  mat W,
                  List S_x,
                  mat sum_nj_xjbar_xjbar_T)
{
    
    int n = size(train_dat)[0];  //データのサンプル数
    int p = size(train_dat)[1];  //ベクトルの次元数
    CharacterVector train_lab_levels = train_lab.attr("levels");  //ラベル
    int G = train_lab_levels.length();  //群の数
    
    // n_g[g,1]でg群の要素数を表す
    int n_g_i = 0;
    for(int i = 0; i < n; ++i){
        if(train_lab(i) == test_lab(0)) n_g_i += 1;
    }
    
    unsigned int i_test_lab = test_lab(0)-1;//要素番号を合わせるためにマイナス１する
    //  NumericMatrix xjxjT_mat0 = xjxjT[i_test_lab];
    mat xjxjT_mat = Rcpp::as<arma::mat>(xjxjT(i_test_lab));
    mat S_x_mat = Rcpp::as<arma::mat>(S_x(i_test_lab));
    
    //CVの外で計算しておいたWやBを使う（群間や群内変動のについて）
    mat bar_x_loo = bar_x;
    //colnames(bar_x_loo) = levels(train_lab)
    bar_x_loo.col(i_test_lab) =  (1 / (double) n_g_i) * ((n_g_i + 1) * bar_x.col(i_test_lab) - test_dat.t());
    mat barx_loo = (1 / (double) n) * ((n + 1) * barx.t() - test_dat);
    
    mat S_x_test1 = xjxjT_mat - test_dat.t() * test_dat;
    mat xjbar_xjbar_T_loo = bar_x_loo.col(i_test_lab) * (bar_x_loo.col(i_test_lab)).t();
    mat S_x_test2 = S_x_test1 - n_g_i * xjbar_xjbar_T_loo;
    mat S_x_test = S_x_test2 / (double) (n_g_i - 1);
    mat W_loo = W - n_g_i * S_x_mat + (n_g_i - 1) * S_x_test;
    
    mat B0 = sum_nj_xjbar_xjbar_T - (n_g_i + 1) * bar_x.col(i_test_lab) * (bar_x.col(i_test_lab)).t() + n_g_i * xjbar_xjbar_T_loo;
    mat B1 = n * barx_loo.t() * barx_loo;
    mat B_loo = B0 - B1;
    
    
    // 固有値・固有ベクトル
    vec eigval;
    mat eigvec;
    mat diagp = eye(p,p);
    eig_sym(eigval, eigvec, 0.5*(W_loo+W_loo.t()) + r * diagp);  // 固有値は小さい順であることに注意
    //eig_sym(eigval, eigvec, W_loo + r * diagp);  // 固有値は小さい順であることに注意、上記は数値誤差を解消するため
    mat Winv_half = eigvec * diagmat(pow(abs(eigval), -0.5)) * eigvec.t();
    eig_sym(eigval, eigvec, 0.5* Winv_half * (B_loo+B_loo.t()) * Winv_half.t());
    //eig_sym(eigval, eigvec, Winv_half * B_loo * Winv_half);
    mat A = Winv_half * eigvec; //射影ベクトル A[,d]でd番目の固有ベクトル
    vec Val = eigval;
    
    //射影先plot(t(Z))で圧縮してプロットできる
    mat Z = (A.cols(p-N_press,p-1)).t() * (test_dat - barx_loo).t(); // 固有値は小さい順なので
    mat barx_loo_rep(p, G, fill::zeros);
    for(int i = 0; i < G; ++i){
        barx_loo_rep.col(i) = barx_loo.t();
    }
    mat M = (A.cols(p-N_press,p-1)).t() * (bar_x_loo - barx_loo_rep);
    mat Z_rep(N_press,G,fill::zeros);
    for(int i = 0; i < G; ++i){
        Z_rep.col(i) = Z;
    }
    vec ones_Npress(N_press, fill::ones);
    vec d = pow(Z_rep - M, 2).t() * ones_Npress;
    CharacterVector hat_lab = Rcpp::as<Rcpp::CharacterVector>(train_lab_levels(d.index_min()));
    bool err = i_test_lab != d.index_min();
    return Rcpp::List::create(Rcpp::Named("err")=err,
                              Rcpp::Named("hat_lab")=hat_lab);   
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List fast_loo_CLD_Rcpp(
                       int p,
                       int cv,
                       IntegerVector test_lab,
                       mat test_dat,
                       IntegerVector train_lab,
                       double hii,
                       rowvec hat_Y_rowvec,
                       rowvec Y_rowvec,
                       mat alpha,
                       vec ci_colvec,
                       vec sum_hatY2,
                       rowvec sum_hY_H_rowvec,
                       double sum_H2_double,
                       mat sum_hatY_j,
                       vec sum_H_j_i,
                       mat bar_x,
                       int N_press,
                       double r
                       )
{
    
    int n = train_lab.length(); //サンプルサイズ
    CharacterVector train_lab_levels = train_lab.attr("levels");  //ラベル
    int G = train_lab_levels.length();  //群の数
    int D = std::min(p, G-1);
    int N_press_adj = std::min(D, N_press);
    uvec n_g(G, fill::zeros);
    for(int i = 0; i < n; ++i){
        n_g(train_lab(i)-1) += 1;
    }
    int n_g_i = n_g(test_lab(0)-1);
    
    mat bar_x_loo = bar_x;
    bar_x_loo.col(test_lab(0)-1) = 1 / (double) n_g_i * ((n_g_i + 1) * bar_x.col(test_lab(0)-1) - vectorise(test_dat));
    
    
    //double hii = diagH(cv);
    rowvec yi_yi_1_h = (hat_Y_rowvec - Y_rowvec)/(1.0 - hii);
    mat alpha_star = alpha + ci_colvec * yi_yi_1_h;
    rowvec sum_x_alpha = (sum_hatY2.t() - pow(hat_Y_rowvec,2)) +
    2.0 * yi_yi_1_h % (sum_hY_H_rowvec - hii * hat_Y_rowvec) +
    (sum_H2_double - hii*hii) * pow(yi_yi_1_h,2);
    vec sum_x_alpha_n = 1.0 / (double) n * sum_x_alpha.t();
    mat beta_star = alpha_star.rows(1,p);
    vec betaTbeta_star = diagvec(beta_star.t() * beta_star);
    vec lambda_star_plus1 = 1.0 / ( sum_x_alpha_n + r * betaTbeta_star );
    
    rowvec tilxi_alpha = 1/(1 - hii) * (hat_Y_rowvec - hii*Y_rowvec);
    mat sum_hatY_j_i = sum_hatY_j;
    sum_hatY_j_i.row(test_lab(0)-1) = sum_hatY_j.row(test_lab(0)-1) - hat_Y_rowvec;
    //vec sum_H_j_i = sum_H_j.col(cv);
    sum_H_j_i(test_lab(0)-1) = sum_H_j_i(test_lab(0)-1) - hii;
    mat ngmat(G, D, fill::zeros);
    for(int i = 0; i < G; ++i){
        for(int j = 0; j < D; ++j){
            ngmat(i,j) = 1.0 / n_g(i);
        }
    }
    mat bar_x_alpha = ngmat % (sum_hatY_j_i + sum_H_j_i * yi_yi_1_h);
    
    
    mat onesmat_G(1,G,fill::ones);
    mat til_bar_x = join_cols(onesmat_G, bar_x_loo);
    mat tilxi_alpha_mat(G,D,fill::zeros);
    mat lambda_star_plus1_mat(G,D,fill::zeros);
    for(int i = 0; i < G; ++i){
        for(int j = 0; j < D; ++j){
            tilxi_alpha_mat(i,j) = tilxi_alpha(j);
            lambda_star_plus1_mat(i,j) = lambda_star_plus1(j)*lambda_star_plus1(j); //廣瀬変更
        }
    }
    mat xi_barx_alpha2 = pow(bar_x_alpha - tilxi_alpha_mat, 2);
    mat lam_xi = lambda_star_plus1_mat % xi_barx_alpha2;
    mat d = lam_xi.cols(0,N_press_adj-1);
    vec ones_Npress(N_press_adj, fill::ones);
    vec d_onesNpress = d * ones_Npress;
    int hat_test_lab_index = d_onesNpress.index_min();
    //# 誤判別率 間違いであれば１増える
    bool err = hat_test_lab_index != test_lab(0)-1;
    
    CharacterVector hat_lab = Rcpp::as<Rcpp::CharacterVector>(train_lab_levels(hat_test_lab_index));
    
    return Rcpp::List::create(Rcpp::Named("err")=err,
                              Rcpp::Named("hat_lab")=hat_lab);
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List cv_Rcpp(mat datx,
             IntegerVector lab,
             double r,
             int N_press,
             List SS_prime,
             CharacterVector type)
{
    
    int N = size(datx)[0];
    int p = size(datx)[1];
    
    mat train_dat = mat(N-1,p,fill::zeros);
    mat test_dat = mat(1,p,fill::zeros);
    IntegerVector train_lab(N-1);
    IntegerVector test_lab(1);
    //CharacterVector train_lab_char_levels(N-1);
    //CharacterVector test_lab_char_levels(1);
    train_lab.attr("class") = "factor";
    test_lab.attr("class") = "factor";
    train_lab.attr("levels") = lab.attr("levels");
    test_lab.attr("levels") = lab.attr("levels");
    
    uvec cv_minus_index = uvec(N-1,fill::zeros);
    
    
    //変数追加
    vec diagH;
    mat hat_Y;
    mat Y;
    mat alpha;
    mat ci;
    vec sum_hatY2;
    mat sum_hY_H;
    vec sum_H2;
    mat sum_hatY_j;
    mat sum_H_j;
    mat bar_x;
    vec barx;
    List xjxjT;
    mat W;
    List S_x;
    mat sum_nj_xjbar_xjbar_T;
    //変数追加
    
    
    
    //変数代入
    if(type(0)=="fast"){         
        diagH = Rcpp::as<arma::vec>(SS_prime["diagH"]);
        hat_Y = Rcpp::as<arma::mat>(SS_prime["hat_Y"]);
        Y = Rcpp::as<arma::mat>(SS_prime["Y"]);
        alpha = Rcpp::as<arma::mat>(SS_prime["alpha"]);
        ci = Rcpp::as<arma::mat>(SS_prime["ci"]);
        sum_hatY2 = Rcpp::as<arma::vec>(SS_prime["sum_hatY2"]);
        sum_hY_H = Rcpp::as<arma::mat>(SS_prime["sum_hY_H"]);
        sum_H2 = Rcpp::as<arma::vec>(SS_prime["sum_H2"]);
        sum_hatY_j = Rcpp::as<arma::mat>(SS_prime["sum_hatY_j"]);
        sum_H_j = Rcpp::as<arma::mat>(SS_prime["sum_H_j"]);
        bar_x = Rcpp::as<arma::mat>(SS_prime["bar_x"]);
    }
    if(type(0)=="exact"){
        bar_x = Rcpp::as<arma::mat>(SS_prime["bar.x"]);
        barx = Rcpp::as<arma::vec>(SS_prime["barx"]);
        xjxjT = Rcpp::as<Rcpp::List>(SS_prime["xjxjT"]);
        W = Rcpp::as<arma::mat>(SS_prime["W"]);
        S_x = Rcpp::as<Rcpp::List>(SS_prime["S.x"]);
        sum_nj_xjbar_xjbar_T = Rcpp::as<arma::mat>(SS_prime["sum_nj_xjbar_xjbar_T"]);
    }
    //変数代入
    
    
    Rcpp::List result;
    double err=0.0;
    for(int cv = 0; cv < N; ++cv){
        int tmp_index=0;
        for(int i = 0; i < N; ++i){
            if(i != cv){
                cv_minus_index(tmp_index) = i;
                train_lab(tmp_index) = lab(i);
                //train_lab_char_levels(tmp_index) = lab_levels(i);
                tmp_index += 1;
            }
        }
        test_dat = datx.row(cv);
        train_dat = datx.rows(cv_minus_index);
        test_lab(0) = lab(cv);
        //train.lab = lab.for[-cv]
        
        //#一個抜きペアリング済みデータで正準判別
        if(type(0)=="exact"){
            
            result = loo_CLD_Rcpp(test_dat,
                                  test_lab,
                                  train_dat,
                                  train_lab,
                                  r,
                                  N_press,
                                  bar_x,
                                  barx,
                                  xjxjT,
                                  W,
                                  S_x,
                                  sum_nj_xjbar_xjbar_T
                                  );
            
        }
        
        if(type(0)=="fast"){
            
            
            result = fast_loo_CLD_Rcpp(p,
                                       cv,
                                       test_lab,
                                       test_dat,
                                       train_lab,
                                       diagH(cv),
                                       hat_Y.row(cv),
                                       Y.row(cv),
                                       alpha,
                                       ci.col(cv),
                                       sum_hatY2,
                                       sum_hY_H.row(cv),
                                       sum_H2(cv),
                                       sum_hatY_j,
                                       sum_H_j.col(cv),
                                       bar_x,
                                       N_press,
                                       r);
            
            
        }
        if(Rcpp::as<Rcpp::LogicalVector>(result["err"])(0)==true) {
            err += 1.0;
        } //#誤判別数を計上
    }
    return Rcpp::List::create(Rcpp::Named("err")=err);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List cv_clustering_Rcpp(mat datx,
                        IntegerVector lab,
                        IntegerVector lab_for,
                        double r,
                        int N_press,
                        List SS_for,
                        CharacterVector Cluster,
                        CharacterVector type
                        )
{
    
    int N = size(datx)[0];
    int p = size(datx)[1];
    
    
    CharacterVector lab_levels = lab.attr("levels");
    CharacterVector lab_for_levels = lab_for.attr("levels");
    
    mat train_dat = mat(N-1,p,fill::zeros);
    mat test_dat = mat(1,p,fill::zeros);
    IntegerVector train_lab(N-1);
    IntegerVector test_lab(1);
    //    CharacterVector train_lab_char_levels(N-1);
    //    CharacterVector test_lab_char_levels(1);
    train_lab.attr("class") = "factor";
    test_lab.attr("class") = "factor";
    train_lab.attr("levels") = lab_for_levels; //どのラベルで定義する注意しておく
    test_lab.attr("levels") = lab_for_levels; //どのラベルで定義する注意しておく
    
    
    int Gall = Cluster.length(); //すべての群の数
    IntegerVector test_lab_clu(1);
    test_lab_clu.attr("class") = "factor";
    uvec lab_uvec = Rcpp::as<arma::uvec>(lab);
    
    //CharacterVector train_lab_char(N-1);  //ラベル
    //CharacterVector train_lab_levels = train_lab.attr("levels");  //ラベル
    
    mat train_dat_clu;
    mat dat_clu;
    IntegerVector lab_clu_export;
    uvec cv_minus_index = uvec(N-1,fill::zeros);
    
    Rcpp::Environment base("package:base");
    Rcpp::Function vector = base["vector"];
    Rcpp::Function is_element = base["is.element"];
    Rcpp::List result;
    Rcpp::List result_clu;
    Rcpp::List SS_for_cv = vector(Rcpp::Named("length")=Gall, Rcpp::Named("mode")="list");
    SS_for_cv.names() = Cluster;
    
    
    //変数追加
    vec diagH, diagH_clu;
    mat hat_Y, hat_Y_clu;
    mat Y, Y_clu;
    mat alpha, alpha_clu;
    mat ci, ci_clu;
    vec sum_hatY2, sum_hatY2_clu;
    mat sum_hY_H, sum_hY_H_clu;
    vec sum_H2, sum_H2_clu;
    mat sum_hatY_j, sum_hatY_j_clu;
    mat sum_H_j, sum_H_j_clu;
    mat bar_x, bar_x_clu;
    vec barx, barx_clu;
    List xjxjT, xjxjT_clu;
    mat W, W_clu;
    List S_x, S_x_clu;
    mat sum_nj_xjbar_xjbar_T, sum_nj_xjbar_xjbar_T_clu;
    //変数追加
    
    
    
    
    //変数代入
    if(type(0)=="fast"){         
        diagH = Rcpp::as<arma::vec>(SS_for["diagH"]);
        hat_Y = Rcpp::as<arma::mat>(SS_for["hat_Y"]);
        Y = Rcpp::as<arma::mat>(SS_for["Y"]);
        alpha = Rcpp::as<arma::mat>(SS_for["alpha"]);
        ci = Rcpp::as<arma::mat>(SS_for["ci"]);
        sum_hatY2 = Rcpp::as<arma::vec>(SS_for["sum_hatY2"]);
        sum_hY_H = Rcpp::as<arma::mat>(SS_for["sum_hY_H"]);
        sum_H2 = Rcpp::as<arma::vec>(SS_for["sum_H2"]);
        sum_hatY_j = Rcpp::as<arma::mat>(SS_for["sum_hatY_j"]);
        sum_H_j = Rcpp::as<arma::mat>(SS_for["sum_H_j"]);
        bar_x = Rcpp::as<arma::mat>(SS_for["bar_x"]);
    }
    if(type(0)=="exact"){
        bar_x = Rcpp::as<arma::mat>(SS_for["bar.x"]);
        barx = Rcpp::as<arma::vec>(SS_for["barx"]);
        xjxjT = Rcpp::as<Rcpp::List>(SS_for["xjxjT"]);
        W = Rcpp::as<arma::mat>(SS_for["W"]);
        S_x = Rcpp::as<Rcpp::List>(SS_for["S.x"]);
        sum_nj_xjbar_xjbar_T = Rcpp::as<arma::mat>(SS_for["sum_nj_xjbar_xjbar_T"]);
    }
    //変数代入
    
    //変数追加
    int index_lab=-1;
    int index_lab_previous=-2;
    List SS_for_cv_List;
    //変数追加
    
    
    
    double err=0.0;
    for(int cv = 0; cv < N; ++cv){
        int tmp_index=0;
        for(int i = 0; i < N; ++i){
            if(i != cv){
                cv_minus_index(tmp_index) = i;
                train_lab(tmp_index) = lab_for(i);
                //train_lab_char_levels(tmp_index) = lab_levels(i);
                tmp_index += 1;
            }
        }
        test_dat = datx.row(cv);
        train_dat = datx.rows(cv_minus_index);
        test_lab(0) = lab_for(cv);
        //train.lab = lab.for[-cv]
        
        //#一個抜きペアリング済みデータで正準判別
        if(type(0)=="exact"){
            result = loo_CLD_Rcpp(test_dat,
                                  test_lab,
                                  train_dat,
                                  train_lab,
                                  r,
                                  N_press,
                                  bar_x,
                                  barx,
                                  xjxjT,
                                  W,
                                  S_x,
                                  sum_nj_xjbar_xjbar_T
                                  );
        }
        if(type(0)=="fast"){
            //result = fast_loo_CLD2_Rcpp(p,
            //                            cv,
            //                            test_lab,
            //                            test_dat,
            //                            train_lab,
            //                            diagH,
            //                            hat_Y,
            //                            Y,
            //                            alpha,
            //                            ci,
            //                            sum_hatY2,
            //                            sum_hY_H,
            //                            sum_H2,
            //                            sum_hatY_j,
            //                            sum_H_j,
            //                            bar_x,
            //                            N_press,
            //                            r);
            
            result = fast_loo_CLD_Rcpp(p,
                                       cv,
                                       test_lab,
                                       test_dat,
                                       train_lab,
                                       diagH(cv),
                                       hat_Y.row(cv),
                                       Y.row(cv),
                                       alpha,
                                       ci.col(cv),
                                       sum_hatY2,
                                       sum_hY_H.row(cv),
                                       sum_H2(cv),
                                       sum_hatY_j,
                                       sum_H_j.col(cv),
                                       bar_x,
                                       N_press,
                                       r);
        }
        
        //#test.labが今回選んだCluster かつ result$hat.lab==test.lab(result$err==FALSE)のときCluster内で正準判別
        if(Rcpp::as<Rcpp::LogicalVector>(is_element(test_lab, Cluster))(0)==true && Rcpp::as<Rcpp::LogicalVector>(result["err"])(0)==false){
            //clu_num.clear();
            //Rprintf("a");
            
            int N_within_cluster=0;
            for(int i = 0; i < N; ++i){
                if(test_lab(0) == lab_for(i)){
                    N_within_cluster += 1;
                }
            }
            uvec clu_num(N_within_cluster);
            uvec clu_num_fortrain(N_within_cluster-1);
            //train_lab_clu.set_size(N_within_cluster-1);
            int tmp_index=0;
            int tmp_index_fortrain=0;
            for(int i = 0; i < N; ++i){
                if(test_lab(0) == lab_for(i)){
                    clu_num(tmp_index) = i;
                    tmp_index += 1;
                    if(i != cv){
                        clu_num_fortrain(tmp_index_fortrain) = i;
                        tmp_index_fortrain += 1;
                    }
                }
            }
            //クラスター内のグループの名前を更新
            uvec levels_withincluster_uvec = unique(lab_uvec(clu_num_fortrain)); //１から始まることに注意
            int levels_size = size(levels_withincluster_uvec)(0);
            CharacterVector levels_withincluster(levels_size);
            for(int i = 0; i < levels_size; ++i){
                levels_withincluster(i) = lab_levels(levels_withincluster_uvec(i)-1);//levels_withincluster_uvec(i)は１からはじまるから
            }
            
            //clu_num = which(lab.for==test.lab)#clu.num[clu.num != cv]でcvを抜いたClusterのdatxでのインデックス
            //test_lab_clu(0) = lab(cv);
            for(int i = 0; i < levels_size; ++i){
                if(levels_withincluster_uvec(i) == (unsigned int)lab(cv)) test_lab_clu(0) = i+1; // test_lab_cluは１からはじまる
            }
            test_lab_clu.attr("levels") = levels_withincluster;
            
            train_dat_clu = datx.rows(clu_num_fortrain);
            IntegerVector train_lab_clu(N_within_cluster-1);
            for(int j = 0; j < N_within_cluster-1; ++j){
                for(int i = 0; i < levels_size; ++i){
                    if(levels_withincluster_uvec(i) == (unsigned int)lab(clu_num_fortrain(j))) train_lab_clu(j) = i+1; // test_lab_cluは１からはじまる
                }
            }
            train_lab_clu.attr("class") = "factor";
            train_lab_clu.attr("levels") = levels_withincluster;
            
            
            //SS.for.cvを更新
            for(int i = 0; i < Gall; ++i){
                if(Cluster(i)==lab_for_levels(test_lab(0)-1)) index_lab = i;
            }
            if(SS_for_cv(index_lab) == R_NilValue){
                dat_clu = datx.rows(clu_num);
                IntegerVector lab_clu(N_within_cluster);
                for(int j = 0; j < N_within_cluster; ++j){
                    for(int i = 0; i < levels_size; ++i){
                        if(levels_withincluster_uvec(i) == (unsigned int)lab(clu_num(j))) lab_clu(j) = i+1; // test_lab_cluは１からはじまる
                    }
                }
                lab_clu.attr("class") = "factor";
                lab_clu.attr("levels") = levels_withincluster;
                //if(type(0)=="exact"){
                //    List a;
                //    a = creasteSS_Rcpp(lab_clu, dat_clu);
                //}
                if(type(0)=="exact") SS_for_cv(index_lab) = creasteSS_Rcpp(lab_clu, dat_clu);
                if(type(0)=="fast") SS_for_cv(index_lab) = creasteSS_fast_Rcpp(lab_clu, dat_clu, r);
                lab_clu_export=lab_clu;
            }
            
            //cvadj
            int cv_adj=-1;
            for(int j = 0; j < N_within_cluster; ++j){
                if(clu_num(j) == (unsigned int)cv) cv_adj=j;
            }
            
            
            
            
            //変数代入
            if(index_lab_previous != index_lab){
                SS_for_cv_List = SS_for_cv[index_lab];
                if(type(0)=="fast"){         
                    diagH_clu = Rcpp::as<arma::vec>(SS_for_cv_List["diagH"]);
                    hat_Y_clu = Rcpp::as<arma::mat>(SS_for_cv_List["hat_Y"]);
                    Y_clu = Rcpp::as<arma::mat>(SS_for_cv_List["Y"]);
                    alpha_clu = Rcpp::as<arma::mat>(SS_for_cv_List["alpha"]);
                    ci_clu = Rcpp::as<arma::mat>(SS_for_cv_List["ci"]);
                    sum_hatY2_clu = Rcpp::as<arma::vec>(SS_for_cv_List["sum_hatY2"]);
                    sum_hY_H_clu = Rcpp::as<arma::mat>(SS_for_cv_List["sum_hY_H"]);
                    sum_H2_clu = Rcpp::as<arma::vec>(SS_for_cv_List["sum_H2"]);
                    sum_hatY_j_clu = Rcpp::as<arma::mat>(SS_for_cv_List["sum_hatY_j"]);
                    sum_H_j_clu = Rcpp::as<arma::mat>(SS_for_cv_List["sum_H_j"]);
                    bar_x_clu = Rcpp::as<arma::mat>(SS_for_cv_List["bar_x"]);
                }
                if(type(0)=="exact"){
                    bar_x_clu = Rcpp::as<arma::mat>(SS_for_cv_List["bar.x"]);
                    barx_clu = Rcpp::as<arma::vec>(SS_for_cv_List["barx"]);
                    xjxjT_clu = Rcpp::as<Rcpp::List>(SS_for_cv_List["xjxjT"]);
                    W_clu = Rcpp::as<arma::mat>(SS_for_cv_List["W"]);
                    S_x_clu = Rcpp::as<Rcpp::List>(SS_for_cv_List["S.x"]);
                    sum_nj_xjbar_xjbar_T_clu = Rcpp::as<arma::mat>(SS_for_cv_List["sum_nj_xjbar_xjbar_T"]);
                }
            }
            //変数代入            
            
            
            //実行
            if(type(0)=="exact"){
                result_clu = loo_CLD_Rcpp(test_dat,
                                          test_lab_clu,
                                          train_dat_clu,
                                          train_lab_clu,
                                          r,
                                          N_press,
                                          bar_x_clu,
                                          barx_clu,
                                          xjxjT_clu,
                                          W_clu,
                                          S_x_clu,
                                          sum_nj_xjbar_xjbar_T_clu
                                          );
            }
            if(type(0)=="fast"){
                //     result_clu = fast_loo_CLD2_Rcpp(p,
                //                                     cv_adj,
                //                                     test_lab_clu,
                //                                     test_dat,
                //                                     train_lab_clu,
                //                                     diagH_clu,
                //                                     hat_Y_clu,
                //                                     Y_clu,
                //                                     alpha_clu,
                //                                     ci_clu,
                //                                     sum_hatY2_clu,
                //                                     sum_hY_H_clu,
                //                                     sum_H2_clu,
                //                                     sum_hatY_j_clu,
                //                                     sum_H_j_clu,
                //                                     bar_x_clu,
                //                                     N_press,
                //                                     r);
                
                result_clu = fast_loo_CLD_Rcpp(p,
                                               cv_adj,
                                               test_lab_clu,
                                               test_dat,
                                               train_lab_clu,
                                               diagH_clu(cv_adj),
                                               hat_Y_clu.row(cv_adj),
                                               Y_clu.row(cv_adj),
                                               alpha_clu,
                                               ci_clu.col(cv_adj),
                                               sum_hatY2_clu,
                                               sum_hY_H_clu.row(cv_adj),
                                               sum_H2_clu(cv_adj),
                                               sum_hatY_j_clu,
                                               sum_H_j_clu.col(cv_adj),
                                               bar_x_clu,
                                               N_press,
                                               r);
                
                
            }
            if(Rcpp::as<Rcpp::LogicalVector>(result_clu["err"])(0)==true) {
                err += 1.0;
            }
            //#誤判別数を計上
            index_lab_previous = index_lab;//廣瀬追加            
        }else{
            if(Rcpp::as<Rcpp::LogicalVector>(result["err"])(0)==true) {
                err += 1.0;
            } //#誤判別数を計上
            
        }
    }
    return Rcpp::List::create(Rcpp::Named("err")=err);
}
