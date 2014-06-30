#include <Rcpp.h>

using namespace Rcpp ;


// struct Deps {
//   std::vector<int> module;
//   double s;
//   double sh;
// };

// std::multimap<int, Deps> restrictTable;


// FIXME: check later penalty for using a string for typeDep
struct geneDeps {
  int typeDep; // smaller, predictable size. A lot less readable, though.
  double s;
  double sh;
  // std::vector<int> deps;
  std::vector< std::vector<int> > deps;
};

// // [[Rcpp::export]]
// void f4(){
//   std::vector<geneDeps> g3;
//   g3.resize(4);
//   g3[0].s = 0.3;
//   g3[0].sh = 0.6;
//   g3[0].deps.resize(2);
//   // g3[0].deps[0].push_back(std::vector<int>(3));
//   std::cout << "g3[0].deps.size() " << g3[0].deps.size() << std::endl;
//   std::cout << "g3[0].deps[0].size() " << 
//     g3[0].deps[0].size() << std::endl;
//   g3[0].deps[0].push_back(3);
//   std::cout << "g3[0].deps[0].size() " << 
//     g3[0].deps[0].size() << std::endl;

//   std::vector<int>vv(4, 3);
//   g3[1].deps.push_back(vv);
//   std::cout << "g3[1].deps.size() " << g3[1].deps.size() << std::endl;
//   g3[1].deps.push_back(std::vector<int>(9, 5));
//   std::cout << "g3[1].deps.size() " << g3[1].deps.size() << std::endl;
  
//   for (auto c : g3[1].deps[0])
//     std::cout << c << ' ';
  
//   for (auto c : g3[1].deps[1])
//     std::cout << c << ' ';
// }

void printRestrictTable(const std::vector<geneDeps>& restrictTable) {
  Rcpp::Rcout << "\n **********  Restriction table inside C++ *******" << std::endl;
  Rcpp::Rcout << "\t Size = " << restrictTable.size() << std::endl;
  for(size_t i = 0; i != restrictTable.size(); ++i) {
    Rcpp::Rcout <<"\n\t\t Dependent node " << (i + 1) << std::endl;
    Rcpp::Rcout <<"\t\t\t typeDep = " << restrictTable[i].typeDep << " ";
    Rcpp::Rcout <<"\t s = " << restrictTable[i].s << " ";
    Rcpp::Rcout <<"\t sh = " << restrictTable[i].sh << std::endl;
    // here the code for parent modules
    Rcpp::Rcout << "\t\t\t Number of parent modules or genes = " << 
      restrictTable[i].deps.size() << std::endl;
    for(size_t j = 0; j != restrictTable[i].deps.size(); ++j) {
      Rcpp::Rcout << "\t\t\t\t Module " << (j + 1) << ": ";
      for (auto c : restrictTable[i].deps[j])
       	Rcpp::Rcout << c << ' ';
      Rcpp::Rcout << std::endl;
    }
    Rcpp::Rcout << std::endl;
  }
}

void restrictTable_to_cpp(Rcpp::List rt,
			  std::vector<geneDeps>& restrictTable) { 
  int ndeps;
  restrictTable.resize(rt.size());

  Rcpp::List rt_element;
  Rcpp::List parent_list;
  Rcpp::IntegerVector module;

  for(int i = 0; i != rt.size(); ++i) {
    rt_element = rt[i];
    restrictTable[i].typeDep = rt_element["type"];
    restrictTable[i].s = rt_element["s"];
    restrictTable[i].sh = rt_element["sh"];
    parent_list = rt_element["parent"];
    ndeps = parent_list.size();

    for(int j = 0; j != ndeps; ++j) {
      module = as<IntegerVector>(parent_list[j]);
      restrictTable[i].deps.push_back(Rcpp::as< std::vector<int> >(module));
    }
  }
}


// [[Rcpp::export]]
void wrap_test_rt(Rcpp::List rtR) {
  std::vector<geneDeps> restrictTable;
  restrictTable_to_cpp(rtR, restrictTable);
  printRestrictTable(restrictTable);
}




//  to be used directly from R, since cannot export: second argument of
//  type unknown to Rcpp // [[Rcpp::export]]

// [[Rcpp::export]]
void restrictTable_to_cpp0(Rcpp::List rt){
  std::vector<geneDeps> restrictTable;
  int ndeps;
  restrictTable.resize(rt.size());
  Rcpp::List rt_element;
  Rcpp::List parent_list;
  Rcpp::IntegerVector module;

  for(int i = 0; i != rt.size(); ++i) {
    rt_element = rt[i];
    restrictTable[i].typeDep = rt_element["type"];
    restrictTable[i].s = rt_element["s"];
    restrictTable[i].sh = rt_element["sh"];
    parent_list = rt_element["parent"];
    ndeps = parent_list.size();

    for(int j = 0; j != ndeps; ++j) {
      module = as<IntegerVector>(parent_list[j]);
      restrictTable[i].deps.push_back(Rcpp::as< std::vector<int> >(module));
    }
  }

  printRestrictTable(restrictTable);

}







