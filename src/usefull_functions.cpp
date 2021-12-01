#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector set_rank (NumericVector x) 
{ int a = 0;
  std::sort(x.begin(),x.end());
  std::reverse(x.begin(), x.end());
  NumericVector out (x.size());
  int rank = 1;
  double curval = x[0] ;
  int e = 0;
  for (a = 0; a < x.size();a++)
  {
    if (x[a] < curval ) {rank++; out[e++]= rank; curval = x[a]; }
    else {out[e++] = rank;}
  }
  
  return out;
}




// [[Rcpp::export]]
Rcpp::DataFrame From_binary_to_transactions (std::vector<std::vector<short>> Tab, std::vector<std::string> name, char sep)
{ int nb_var = Tab.size();
  int nb_ind = Tab[0].size();
  std::vector<std::string> out;
  out.reserve(nb_ind);
  int j = 0;
  for (int i=0; i < nb_ind;i++)
  { std::string str; 
    for (j=0; j< nb_var;j++)
    {if (Tab[j][i]) 
    {str += name[j];
      str += sep;}}
    str.erase(str.size()-1);
    out.push_back(str);}
  
  return Rcpp::DataFrame::create(Rcpp::Named("Transactions")=out);
  
}

// [[Rcpp::export]]
Rcpp::DataFrame From_transaction_to_binary( std::vector<std::string> transac, char deli)
{ int nb_individus = transac.size();
  
  std::vector <std::string> listename;
  std::unordered_map<std::string,int> maptr;
  for (std::string s :  transac)
  { std::string q ="";
    for (char c :s)
      { if (c!= deli) {q+= c;}
      else {maptr[q]; q.erase(q.begin(),q.end());}
      }
    maptr[q];
  }
  
  int nb_var = maptr.size();
  
  std::unordered_map<std::string,int>::iterator itmap;
  std::unordered_map<std::string,int>::iterator itendmap = maptr.end();
  int a = 0;
  for (itmap = maptr.begin();itmap != itendmap;itmap++)
  {itmap->second = a++;
    listename.push_back(itmap->first);}
  
  std::vector<std::vector<int>> mytab (nb_var, std::vector<int>(nb_individus));
  int i = 0;
  for (std::string s :  transac)
  { std::string q ="";
    for (char c :s)
      { if (c!= deli) {q+= c;}
      else {mytab[maptr[q]][i]=1; q.erase(q.begin(),q.end());}
      }
    mytab[maptr[q]][i++]=1;
  }
  
  Rcpp::DataFrame Dataout (mytab);
  Dataout.attr("names") = listename;
  
  return Dataout;
}


int check_target (std::string &set, std::vector<std::string> & target_item, int size_tgt)
{ 
  int ret = 0; int i,j,k,size_set = 0;
  for (char c : set) {if (c==',') size_set++;}
  k = 0;
  std::string stemp = "";
  std::vector<std::string> unique_item;
  unique_item.reserve(size_set);
  for (char c: set)
  {
    if (c ==',') 
    { stemp +=c;
      unique_item.push_back(stemp);
      
      stemp = ""; continue;
    }
    stemp += c;
  }
  k = 0;

  for (std::string st : unique_item)
  { 
    for (i = 0; i < size_tgt;i++) {if (target_item[i] == st) k++; }
    
  }
  
  if (k >= size_tgt) {return 1;}
  return 0;
}

void set_tgt_vec_unique (std::vector<std::string> &tgt_vec_unique, std::string tgt)
{ std::string stemp;
  int i = 0;
  for (char c: tgt)
  {
    if (c ==',') 
    { stemp +=c;
      tgt_vec_unique[i++] = stemp;
      
      stemp = ""; continue;
    }
    stemp += c;
  }
  
}

// [[Rcpp::export]]
Rcpp::DataFrame is_antecedant_or_consequent ( std::vector<std::string> set_to_check, std::string target, std::string ante_or_cons )
{
  int nb_elem = set_to_check.size();
  std::string type_tg = ante_or_cons ;
  if (ante_or_cons == "a") {type_tg = "antecedant";}
  if (ante_or_cons == "c") {type_tg = "consequent";}
  int size_tg =0;
  for (char c : target) {if (c ==',') size_tg++;}
  std::vector<std::string> unique_item_in_target (size_tg) ;
  set_tgt_vec_unique(unique_item_in_target, target);
  NumericVector out (nb_elem);
  for (int i = 0; i < nb_elem;i++)
  {out[i] = check_target(set_to_check[i],unique_item_in_target,size_tg);}
  std::string name_out = target +"_is_";
  name_out += type_tg;
  DataFrame bool_target = Rcpp::DataFrame::create(Rcpp::Named(name_out)=out);
  return bool_target;
  
}  