#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <fstream>
#include <math.h>

Environment base = Environment("package:base");

struct pnodesr
{
  pnodesr (double sup_, int litem_, int size_, int is_target_ ) : sup(sup_), litem(litem_), size(size_), is_target(is_target_){};
  pnodesr * brother = NULL;
  pnodesr * son = NULL;
  pnodesr * father = NULL; 
  double sup = 0;
  bool is_target = 0;
  int size =0;
  int litem =0;
};
pnodesr * clist = NULL;
pnodesr * ctree = NULL;
pnodesr * cfath = NULL;


struct rules
  { rules (double conf_, std::string  ante_, std::string cons_, std::string set_) : conf(conf_), ante(ante_), cons(cons_), set(set_) {};
    double conf =0;
    std::string ante ;
    std::string cons ;
    std::string set ;
    bool valid = 0;
  };


struct fction_selec
{
  double (*pfunction)(double & freq_ante, double & freq_complem, double & freq_set );
  double treshold = 0;
  double max = 0;
  std::string name_fction = "";
  
};


int st_toi (std::string coeff)
{
  int i = 0;
  for (char c : coeff) {i+= (int)c;}
  return i;
}
  


double conftwo (double & freq_ante, double & freq_complem, double & freq_set )
{
  return freq_set/freq_ante;
  
}

double powertwo (double & freq_ante, double & freq_complem, double & freq_set )
{
  return (1-freq_ante - freq_complem + freq_set)/(1-freq_ante);
}

double confone (double & freq_ante, double & freq_complem, double & freq_set )
{
  return (1-freq_ante - freq_complem + freq_set)/(1-freq_complem);
}

double powerone (double & freq_ante, double & freq_complem, double & freq_set )
{
  return freq_set/freq_complem;
}

double cov (double & freq_ante, double & freq_complem, double & freq_set )
{
  return freq_set - freq_ante*freq_complem;
}

double corr (double & freq_ante, double & freq_complem, double & freq_set)
{ 
  double cov = freq_set - freq_ante*freq_complem;
  return cov/std::pow(freq_ante*(1-freq_ante)*freq_complem*(1-freq_complem),0.5);
}


double kappa (double & freq_ante, double & freq_complem, double & freq_set)
{ 
  double cov = freq_set - freq_ante*freq_complem;
  return cov/(freq_ante*(1-freq_complem)+freq_complem*(1-freq_ante));
}

double maxwellP (double & freq_ante, double & freq_complem, double & freq_set)
{
  double cov = freq_set - freq_ante*freq_complem;
  return (2*cov)/(freq_ante*(1-freq_ante)+freq_complem*(1-freq_complem));
}

void set_coeff_fc (fction_selec * tab_fcs, int * tab_map , std::vector<std::string> & complete_coeff ,int nb_coeff)
{ 
  for (int i = 0; i < nb_coeff;i++)
  { 
    switch (tab_map[i])
    {
    case 0:
      tab_fcs[i].name_fction = complete_coeff[0];
      tab_fcs[i].pfunction = &conftwo;
      break;

    case 1:
      tab_fcs[i].name_fction = complete_coeff[1];
      tab_fcs[i].pfunction = &confone;
      break;

    case 2:
      tab_fcs[i].name_fction = complete_coeff[2];
      tab_fcs[i].pfunction = &powertwo;
      break;
     
      
    case 3:
      tab_fcs[i].name_fction = complete_coeff[3];
      tab_fcs[i].pfunction = &powerone;
      break;

      
    case 4:
      tab_fcs[i].name_fction = complete_coeff[4];
      tab_fcs[i].pfunction = &cov;
      break;

      
    case 5:
      tab_fcs[i].name_fction = complete_coeff[5];
      tab_fcs[i].pfunction = &corr;
      break;

      
    case 6:
      tab_fcs[i].name_fction = complete_coeff[6];
      tab_fcs[i].pfunction = &kappa;
      break;

      
    case 7:
      tab_fcs[i].name_fction = complete_coeff[7];
      tab_fcs[i].pfunction = &maxwellP;
      break;
      
    default :
      tab_fcs[i].name_fction = complete_coeff[0];
      tab_fcs[i].pfunction = &conftwo;      
    }
  }   
  
  
}  


int check_target (std::string &set, std::vector<std::string> &targets)
{
  int ret = 0; 
  std::string stemp = "";
  for (char c: set)
  {
    if (c ==',') 
    { stemp +=c;
      
      for (std::string & st : targets)
      { 
        if (stemp == st) ret =1;
      }
      stemp = ""; continue;
    }
    stemp += c;
  }

  return ret;
}


void set_fction_select (fction_selec * tab_fcs, fction_selec * tab_fcs_choice ,std::vector<std::string> & coeff, std::vector<std::string> & coeff_choice ,std::vector<double> & val_coeff)
  { std::vector<std::string> complete_coeff = {"Conf2","Conf1","Power2","Power1","Cov","Corr","Kappa","Maxp"};
    std::unordered_map<std::string,int> map_abuse;
    int i;
    int nb_coeff = coeff.size();
    int nb_coeff_choice = coeff_choice.size();
    for (i = 0; i < complete_coeff.size();i++) {map_abuse[complete_coeff[i]] = i; }

    int tab_sok [nb_coeff];
    int tab_schoice [nb_coeff_choice];
    for ( i = 0; i < nb_coeff;i++)
    {
      tab_sok [i]= map_abuse[coeff[i]];
      tab_fcs [i].treshold = val_coeff[i];

    }
    
    for ( i = 0; i < nb_coeff_choice;i++)
    {
      tab_schoice [i]= map_abuse[coeff_choice[i]];
    }
    set_coeff_fc(tab_fcs,tab_sok,complete_coeff,nb_coeff);
    set_coeff_fc(tab_fcs_choice,tab_schoice,complete_coeff,nb_coeff_choice);
  }



void set_target_indices ( bool * ind_target, std::vector<int> & litem_glob , std::vector<int> & size_glob, std::vector<int> & litem_target, int & nb_it)
{ 
  int i,j,k;
  Rcpp::Rcout << "nb_it = " << nb_it << "\n";
  for(i = 0; i < nb_it;i++){ind_target[i]=0;}
  nb_it = 0;
  int v = 1;
  int nb_g = litem_glob.size();
  int nb_t = litem_target.size();
  litem_glob.push_back(-1);
  int cursize = 1;
  int curlitem = 0;
  int tab_h [nb_g];
  for (i = 0; i < nb_g;)
  {curlitem = litem_glob[i];
    v=1;
    for (j =0; j <  nb_t;j++)
    { 
      if (litem_target[j]==curlitem)
      { cursize = size_glob[i];
       
        tab_h[nb_it++] = i;
        v=0;
        i++;
        while (size_glob[i] > cursize)
        {
          tab_h[nb_it++] = i;
          i++;
        }
        
        break;
      }
    }
    if(v!=0) i++;
  }
  for ( i =0; i < nb_it;i++) {ind_target[tab_h[i]] = 1;}
}

void set_target_labels (std::vector<std::string> & nameliste, std::vector<std::string> & targets, std::vector<int> & labeliste, std::vector<int> & label_targets)
{ int r,t;
  int name_sl = nameliste.size();
  int nb_target = targets.size();
  
  for (r = 0; r < nb_target;r++)
  { std::string tgt = targets[r];
    Rcpp::Rcout << "target est " << tgt << "\n";
    for ( t = 0;t < name_sl;t++ )
    {
      if (nameliste[t] == tgt )  {label_targets.push_back(t); Rcpp::Rcout << "le label est "<< t << "\n"; break;}
    }
  }
}

int give_size_set (std::string & str )
{ int s = 0;
  for (char c : str) {if (c ==',') s++; }
  return s;
}

int check_rules (fction_selec * tab_fction, int nb_fction, std::list<rules> & ruleliste, std::unordered_map<std::string,double> & map_sup)
{  
    
    int i,valid;
    int nb_valid = 0;
    double freq_ante,freq_cons,freq_set;
    std::list<rules>::iterator ite = ruleliste.end();
    std::list<rules>::iterator itb = ruleliste.begin();

     
   for (itb = ruleliste.begin(); itb!=ite;itb++)
   { valid = 1;
     
     freq_ante = map_sup[itb->ante]; freq_cons = map_sup[itb->cons]; freq_set = map_sup[itb->set];
     //Rcpp::Rcout << freq_ante << " " << freq_cons << " " << freq_set << "\n";
      for (i = 0; i < nb_fction;i++)
      { //Rcpp::Rcout << conftwo(freq_ante,freq_cons,freq_set) << "\n";   
          if (tab_fction[i].pfunction(freq_ante,freq_cons,freq_set)<tab_fction[i].treshold){valid = 0; break;}
        }
      if (valid) {itb->valid=1; nb_valid++;}
    }
    return nb_valid;
} 

void fill_rules (fction_selec * tab_fction, int nb_fction, std::list<rules> & ruleliste, std::unordered_map<std::string,double> & map_sup, std::vector<std::vector<double>>&  Mat_coeff,std::vector<std::vector<int>>& Mat_size, std::vector<std::vector<std::string>> & Mat_names, bool ok_sup, bool ok_size )
{ int i;int sit = 0;

  double freq_ante,freq_cons,freq_set;
  std::list<rules>::iterator ite = ruleliste.end();
  std::list<rules>::iterator itb = ruleliste.begin();

  for (itb = ruleliste.begin(); itb!=ite;itb++)
  { if(!itb->valid) continue;
    freq_ante = map_sup[itb->ante]; freq_cons = map_sup[itb->cons]; freq_set = map_sup[itb->set];
    Mat_names[0][sit] = itb->ante; Mat_names[1][sit] = itb->cons; Mat_names[2][sit] = itb->set;
    if(ok_size) {Mat_size[0][sit] = give_size_set(itb->ante); Mat_size[1][sit] = give_size_set(itb->cons); Mat_size[2][sit] = give_size_set(itb->set);}
    for(i = 0; i < nb_fction;i++)
      {
        Mat_coeff[i][sit] = tab_fction[i].pfunction(freq_ante,freq_cons,freq_set);
      }
    if(ok_sup) {Mat_coeff[nb_fction][sit] = freq_ante; Mat_coeff[nb_fction+1][sit]=freq_cons; Mat_coeff[nb_fction+2][sit]=freq_set;}
    sit++;  
   
  }
  
} 


void gen_tree (pnodesr * curlist, pnodesr * curtree, pnodesr * curfather,pnodesr ** listepnodesr, int  & sit, int & stop)
{ if (sit != stop)
{
  if (curlist->size > curtree->size)
  { curlist->father = curtree; 
    curtree->son = curlist;
    sit++;
    gen_tree (listepnodesr[sit],curlist,curtree,listepnodesr,sit,stop);
  }
  if (curlist->size == curtree->size)
  {   curlist->father = curfather;
    curtree->brother = curlist;
    sit++;
    gen_tree (listepnodesr[sit],curlist,curfather,listepnodesr,sit,stop);
    
  }
  
  if (curlist->size < curtree->size)
  {
    gen_tree (curlist, curtree->father,curtree->father->father,listepnodesr,sit,stop);
  }
}
else {clist = curlist; ctree = curtree; cfath = curfather;}   

}


void rebuild_tree (std::vector<double> & supvalue, std::vector<int> & sizevalue, std::vector<int> & litemvalue, bool * target_tab ,pnodesr *& rootrules, pnodesr** tabpnodes ,int nb_freq )
{
  int sit = 0; int nbg = 0; int j; 
  int k = nb_freq/10000;
  int reste_k = nb_freq-10000*k;
  int ct = 1;


  for (j = 0; j < nb_freq;j++)
    { 
      tabpnodes[j] = new pnodesr(supvalue[j],litemvalue[j],sizevalue[j],target_tab[j]);
    }
  
  clist = tabpnodes[0];
  ctree = rootrules;
  cfath = rootrules;
  for(j = 0; j< k; j++,ct++)
    {
      nbg = 10000 * ct;
      gen_tree(clist,ctree,cfath,tabpnodes,sit,nbg);
    }
  gen_tree(clist,ctree,cfath,tabpnodes,sit,nb_freq);
}




void allrules  (std::string & stree, std::string& cntree , std::string * tabconseqalpha, int size , int good, std::unordered_map<std::string,double>& map_double, double cursup, int cursize, double Minconf,std::list<rules>& ruliste, std::string & strset )
{ 
  
  std::string trstree = stree;
  std::string cnstree = cntree;
  int ngr = 0;
  int tail = size - good;
  
  if (tail == 1)
  {  
    trstree.erase(trstree.find(tabconseqalpha[good]),tabconseqalpha[good].size());
    double indic = conftwo(map_double[trstree],map_double[cnstree],cursup);
    if (indic > Minconf)
    { 
      cnstree +=  tabconseqalpha[good];
      ruliste.emplace_front (rules(indic,trstree,cnstree,strset));
    }
    
    
  }
  else
  { 
    std::string * antecedant = new std::string [tail];
    std::string * consequent  = new std::string [tail];
    std::string * consalpha   = new std::string [tail];
    
    for (int i = good; i < size ; i++)
    {
      
      trstree.erase(trstree.find(tabconseqalpha[i]),tabconseqalpha[i].size());
      double indic = conftwo(map_double[trstree],map_double[cnstree],cursup);
      if (indic > Minconf)
      { 
        cnstree +=  tabconseqalpha[i];
        antecedant[ngr] = trstree;
        consequent[ngr] = cnstree;
        consalpha[ngr++]= tabconseqalpha[i];
        ruliste.emplace_front (rules (indic,trstree,cnstree,strset));
      }
      trstree = stree; 
      cnstree = cntree;
    }
    
    if (ngr >1 && cursize> 2) 
    {  
      int ngt = ngr -1;
      for (int j = 0; j < ngt;j++)
      { 
        allrules(antecedant[j],consequent[j],consalpha,ngr,j+1,map_double,cursup,cursize-1,Minconf, ruliste,strset);
      }
    }
    
    delete [] antecedant;
    delete [] consequent;
    delete [] consalpha;
    
  }
  
  
}


void genrules (pnodesr & tree ,int size ,std::string ** tabname,std::string & str, std::string * varnames, std::unordered_map<std::string,double>& map_double, double& Minconf, std::list<rules>& ruliste)
{ 
  std::string strtemp = str + varnames[tree.litem];
  double cursup = tree.sup;
  map_double[strtemp] = tree.sup;
  tabname[size-1] = &varnames[tree.litem];
  
  if(tree.is_target){
  
  std::string * antecedant = new std::string [size];
  std::string * consequent  = new std::string [size];
  std::string * consalpha   =  new std::string [size];
  int nrg = 0;     
  std::string tmpstr ;
  tmpstr.reserve(100);
  tmpstr +=  *(tabname[1]);
  for (int i = 2; i < size;i++)
  { tmpstr += *(tabname[i]);}
  double indic = conftwo(map_double[tmpstr],map_double[*tabname[0]],cursup);
  if (indic > Minconf ) {
    ruliste.emplace_front(rules(indic,tmpstr,*tabname[0],strtemp));
    antecedant[nrg] = tmpstr;
    consalpha[nrg] = *tabname[0];
    consequent [nrg++] = *tabname[0];}
  
  tmpstr.erase(tmpstr.begin(),tmpstr.end());
  tmpstr += *(tabname[0]);
  
  size_t st = tabname[0]->size();
  int l = 0;
  for (int j = 1; j < size; j++)
  {
    for (l =1; l < size ; l++)
    {
      if (l != j) tmpstr += *(tabname[l]);     
    }
    double indic = conftwo(map_double[tmpstr], map_double[*tabname[j]],cursup);
    if (indic > Minconf)  
    {
      ruliste.emplace_front(rules(indic,tmpstr,*tabname[j],strtemp));
      antecedant[nrg] = tmpstr;
      consalpha [nrg] = *tabname[j];
      consequent[nrg++] = *tabname[j];
    } 
    tmpstr.erase(tmpstr.begin()+st,tmpstr.end());       
  }
  
  if (nrg >1 && size > 2) 
  { int ngt = nrg -1;
    for (int m = 0; m < ngt;m++)
    { 
      allrules(antecedant[m],consequent[m],consalpha,nrg,m+1,map_double,cursup,size-1,Minconf, ruliste,strtemp);
    }
  }
  delete [] antecedant;
  delete [] consequent;
  delete [] consalpha; 
  }
  
  if (tree.son) 
    {
     genrules(*tree.son,size+1,tabname,strtemp,varnames,map_double,Minconf,ruliste);
     }
  
  if (tree.brother)
    { 
      genrules(*tree.brother,size,tabname,str,varnames,map_double,Minconf,ruliste);
    }

} 


void Genrules_root (pnodesr & roots, int size, std::string ** tabname, std::string &str, std::string * varnames, std::unordered_map<std::string,double>& map_double,double & Minconf, std::list<rules>& ruliste)
{ std::string strtemp = str + varnames[roots.litem];
  double cursup = roots.sup;
  map_double[strtemp] = cursup;
  tabname[0] = &varnames[roots.litem];
  if (roots.son)
  { genrules (*roots.son,size+1,tabname,strtemp,varnames,map_double, Minconf,ruliste);}
  if (roots.brother) {  Genrules_root(*roots.brother,size,tabname,str,varnames,map_double, Minconf,ruliste);}    
}



void erase_tree (pnodesr *& Tree)
{ 
  if (Tree->son) {erase_tree(Tree->son);}
  if (Tree->brother) {erase_tree(Tree->brother);}
  delete Tree;
}




// [[Rcpp::export]]
Rcpp::DataFrame prefrulestrat (Rcpp::List indic ,std::vector<std::string> & nameliste ,std::vector<std::string> coeffs, std::vector<double> valcoef, std::vector<std::string> coeff_extract, std::vector<std::string> targets, std::vector<int> ok_choice )
{ Rcpp::Function readline = base["readline"];
  if (coeffs[0]!= "Conf2" || valcoef.size()==0) 
    { Rcpp::Rcout << "Please set Conf2 as first coefficient and set in valcoef its treashold \n ";
      std::vector<int> outnul(1);  
      return Rcpp::DataFrame::create( Rcpp::Named("No Conf2 value")= outnul);}
  
  int choice_quit = 0;
       
  
  double minConf = valcoef[0];
  Rcpp::Rcout << "La valeur de Minconf est " << minConf << "\n";
  bool ok_size = ok_choice[0];
  bool ok_sup = ok_choice[1];
  bool ok_global_indic = ok_choice[2];
  bool ok_targets = 1;
  if (coeff_extract[0] == "all_indicators") {coeff_extract = std::vector<std::string> {"Conf2","Conf1","Power2","Power1","Cov","Corr","Kappa","Maxp"};}
  if (targets[0]=="all_tgts") {ok_targets =0;}
  
  int r = 0; int t = 0;

  int nb_coeff = coeffs.size();
  int nb_coeff_extract = coeff_extract.size();
  std::vector<int> litem_targets;
  int nb_target = targets.size();
  std::vector<double> supvalue =  indic[0]; 
  std::vector<int> sizevalue = indic[1];
  std::vector<int> litemvalue = indic [2];
  int nb_freq = supvalue.size();
  bool * ind_target = new bool [nb_freq];
  for (r = 0; r < nb_freq;r++) {ind_target[r] = 1;}
  Rcpp::Rcout << "nb_freq = " << nb_freq << "\n";
  int nb_target_frequent = nb_freq;
  pnodesr ** tabpnodes = new pnodesr* [nb_freq];
  pnodesr * rootrules =  new pnodesr(0,0,0,0);
  if (ok_targets) 
  {  
  set_target_labels(nameliste,targets,litemvalue,litem_targets);
  set_target_indices(ind_target,litemvalue,sizevalue,litem_targets,nb_target_frequent);
  
  if (nb_freq == 0)
  {
    Rcpp::Rcout << "No one rule with this target(s) \n";
    std::vector<int> outnul(1);  
    return Rcpp::DataFrame::create( Rcpp::Named("No rules target")= outnul);
  }
  }
  rebuild_tree(supvalue,sizevalue,litemvalue,ind_target,rootrules,tabpnodes,nb_freq);
  
  delete [] ind_target;
  
  nb_freq = nb_target_frequent;

  
  std::vector<std::string> complete_coeff = {"Conf2","Conf1","Power2","Power1","Cov","Corr","Kappa","Maxp"};
  
  std::unordered_map<std::string,double> map_sup;
  map_sup.reserve(nb_freq);
  fction_selec tab_choice_coeff [nb_coeff+1];
  tab_choice_coeff[0].pfunction = &conftwo;
  fction_selec tab_choice_extract [nb_coeff_extract+1];
  
  set_fction_select (tab_choice_coeff,tab_choice_extract,coeffs,coeff_extract,valcoef);

  
  
  std::vector<std::string> names_size {"Size_Ante","Size_Cons","Size_Set"};
  std::vector<std::string> names_sup {"Sup_Ante","Sup_Cons","Sup_Set"};

  std::list<rules> ruleliste;
  int taille_nameliste = nameliste.size();
  std::string * tabname  [taille_nameliste];
  std::string st;
  st.erase(st.begin(),st.end());
  st.reserve (100);
  int compte = 0;
  Genrules_root (*rootrules->son,1,tabname,st,&nameliste[0],map_sup, minConf ,ruleliste);
  erase_tree(rootrules);
  

  int n_rules = check_rules(tab_choice_coeff,nb_coeff,ruleliste,map_sup);
  Rcpp::Rcout << "n_rules = " << n_rules << "\n"; 

  std::vector<std::vector<std::string>> Mat_names (3, std::vector<std::string>(n_rules)); 
  std::vector<std::vector<int>> Mat_size;
  std::vector<std::vector<double>> Mat_coeff;
  
  if (ok_size) {Mat_size = std::vector<std::vector<int>> (3, std::vector<int> (n_rules));}
  if (ok_sup && ok_global_indic )  {Mat_coeff = std::vector<std::vector<double>> (nb_coeff_extract+5, std::vector<double> (n_rules));}
  if (ok_sup && !ok_global_indic )  {Mat_coeff = std::vector<std::vector<double>> (nb_coeff_extract+3, std::vector<double> (n_rules));}
  if (!ok_sup && ok_global_indic )  {Mat_coeff = std::vector<std::vector<double>> (nb_coeff_extract+2, std::vector<double> (n_rules));}
  if (!ok_sup && !ok_global_indic )  {Mat_coeff = std::vector<std::vector<double>> (nb_coeff_extract, std::vector<double> (n_rules));}
  

  fill_rules(tab_choice_extract,nb_coeff_extract,ruleliste, map_sup,Mat_coeff,Mat_size,Mat_names,ok_size,ok_sup);

  int nb_frequent_target = 0;
  for (std::unordered_map<std::string,double>::iterator it = map_sup.begin();it != map_sup.end();it++)
  { std::string sttest = it->first;
    nb_frequent_target += check_target(sttest,targets);
  }
  // Rcpp::Rcout << "IL y a " << nb_frequent_target << " frequents avec la cible \n";
  if (ok_global_indic)
    { std::list<rules>::iterator it = ruleliste.begin();
      double * pt_meanc = &Mat_coeff[Mat_coeff.size()-2][0];
      double * pt_strong = &Mat_coeff[Mat_coeff.size()-1][0];
      for (r = 0; r < n_rules; it++)
      { if (!it->valid) continue;
          pt_meanc[r] = (conftwo(map_sup[it->ante],map_sup[it->cons],map_sup[it->set])+confone(map_sup[it->ante],map_sup[it->cons],map_sup[it->set])+powertwo(map_sup[it->ante],map_sup[it->cons],map_sup[it->set])+powerone(map_sup[it->ante],map_sup[it->cons],map_sup[it->set]))/(double)4;
          pt_strong[r] = (double)(1-map_sup[it->ante]-map_sup[it->cons]+2*map_sup[it->set]);
          r++;
        }
    }
  
  std::vector<std::string> names_dataframe {"Antecedant","Consequent","set"};

  
  /*
  std::vector<int> target_antecedant;
  std::vector<int> target_consequent;
  
  if (ok_targets)
  { 
    
    target_antecedant= std::vector<int> (n_rules);
    target_consequent= std::vector<int> (n_rules);
    std::string * pt_tante = &Mat_names[0][0];
    std::string * pt_tcons = &Mat_names[1][0];
    
    for (r =0; r < n_rules;r++)
      {
        target_antecedant[r]=(check_target(pt_tante[r],targets));
        target_consequent[r]=(check_target(pt_tcons[r],targets));
      }
    
  }
  */
  
  
  
    
  Rcpp::List prep_out; for (r = 0; r < Mat_names.size();r++) {prep_out.push_back(Mat_names[r]);}
  if (ok_size)
  {
    for (r = 0; r < names_size.size();r++)
    {
      names_dataframe.push_back(names_size[r]);
      prep_out.push_back(Mat_size[r]);
    }
  }
  if (ok_sup)
  { int m = 0;
      for (r = nb_coeff_extract; r < nb_coeff_extract+3; r++)
        {
          names_dataframe.push_back(names_sup[m++]);
          prep_out.push_back(Mat_coeff[r]);
        }
    }

  
  for (r = 0;r<nb_coeff_extract;r++)
    {
      names_dataframe.push_back(coeff_extract[r]);
      prep_out.push_back(Mat_coeff[r]);
    }
  
  if (ok_global_indic)
  {
    names_dataframe.push_back("meanc");
    names_dataframe.push_back("strong");
    prep_out.push_back(Mat_coeff[Mat_coeff.size()-2]);
    prep_out.push_back(Mat_coeff[Mat_coeff.size()-1]);
  }
  
  /*if (ok_targets)
  {
    names_dataframe.push_back("target_A");
    names_dataframe.push_back("target_B");
    prep_out.push_back(target_antecedant);
    prep_out.push_back(target_consequent);
    
  }*/

  Rcpp::DataFrame out_data (prep_out);
  out_data.attr("names")=names_dataframe;
  
  return out_data;
                                                 
}