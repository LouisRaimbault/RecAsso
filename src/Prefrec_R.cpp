#include <Rcpp.h>
using namespace Rcpp;
#pragma GCC target("popcnt")

Environment base = Environment("package:base");

int nbfreq = 0;
static int minsup =0;
int curlitem = 0;
static int nbc = 0;
static uint64_t Ultab [65];

struct pnodes 
{
  pnodes ( int sup_, int litem_, pnodes* brother_ );
  int sup = 0;
  unsigned int litem =0;
  pnodes * brother = NULL;
  pnodes * son = NULL;
};

struct rnodes
  {   rnodes (){}
    rnodes ( int sup_,bool indic_, uint64_t* tab_);
    int sup = 0;
    bool indic =0;
    uint64_t * tab = NULL;
  };


pnodes::pnodes (int sup_, int litem_, pnodes* brother_ ):sup(sup_), litem(litem_), brother(brother_) {};
rnodes::rnodes ( int sup_,bool indic_, uint64_t* tab_):sup(sup_), indic(indic_), tab(tab_){};

void Init_data (std::vector<std::string> & transactions, int & nrows, int & maxul ,int & nvar, uint64_t **& Bitdata, std::vector<std::string> & varnames, char deli ) 
{ 
  Rcpp::Rcout << "import data from file as char array ... \n";
  char * tabchar;
  long longueur = 0;
  long i = 0;
  char c;
  int sit =0;
  int t = 0;
  std::string st ="";
  st.reserve(1000);

  std::unordered_map<std::string,uint64_t*> maptr;
  st = "";
  for (std::string s :transactions)
  { 
    for (char c :s) 
    {
      if (c!= deli) {st+= c;}            
      else {maptr[st]; st="";}
    }
    maptr[st]; st="";
  }
  Rcpp::Rcout << "Done. Preparing BitData store in uint64_t ... ";
  nvar = maptr.size();
  maxul = nrows/64;
  int rst = nrows-64*maxul;
  if (rst != 0) maxul++;
  Bitdata = new uint64_t * [nvar];
  for (t = 0; t < nvar;t++) {Bitdata[t] = new uint64_t [maxul]{0};}
  
  varnames.reserve(nvar+1);
  std::unordered_map<std::string,uint64_t*>::iterator itmap;
  std::unordered_map<std::string,uint64_t*>::iterator itendmap = maptr.end();
  t=0; i = 0; sit = 0;
  for (itmap = maptr.begin();itmap != itendmap;itmap++)
  {itmap->second = Bitdata[t++];
    varnames.push_back(itmap->first+","); 
  }
  t = 0;
  for (std::string s : transactions)
  { st ="";
    for (char c :s)
      { if (c!= deli) {st+= c;}
      else {maptr[st][t]|= 1UL << i; st="";}
      }
    maptr[st][t]|= 1UL << i;
    sit++; 
    t = sit/64;
    i = sit%64;
  }
  
  Rcpp::Rcout << " ok \n";
}


void sortoutdata (uint64_t** Bdata, int * sumvec, int nb_elem, std::vector<std::string> &varnames, char ordre)
{
  int i,j,treme,itreme,temp;
  uint64_t * tempt;
  std::string tempstr;
  if (ordre=='d')
  {  
    for (i = 0; i<nb_elem-1;i++)
    {itreme = i; treme=sumvec[i];
    for (j=i+1;j<nb_elem;j++) {if (sumvec[j]>treme) {treme = sumvec[j];itreme=j;}}
    temp=sumvec[itreme];sumvec[itreme]=sumvec[i];sumvec[i]=temp;
    tempt = Bdata[itreme]; Bdata[itreme] = Bdata[i];Bdata[i]=tempt;
    tempstr = varnames[itreme]; varnames[itreme] = varnames[i];varnames[i]=tempstr;
    }
  }
  if (ordre == 'i')
  {  
    for (i = 0; i<nb_elem-1;i++)
    {itreme = i; treme=sumvec[i];
    for (j=i+1;j<nb_elem;j++) {if (sumvec[j]<treme) {treme = sumvec[j];itreme=j;}}
    temp=sumvec[itreme];sumvec[itreme]=sumvec[i];sumvec[i]=temp;
    tempt = Bdata[itreme]; Bdata[itreme] = Bdata[i];Bdata[i]=tempt;
    tempstr = varnames[itreme]; varnames[itreme] = varnames[i];varnames[i]=tempstr;
    }
  }
}

int init_prefixtree (uint64_t ** Bitdata, uint64_t **& Bdata,std::vector<std::string>& varnames, int support, int nvar, int maxul, char ordre )
{ 
  int sumvec [nvar];
  int tokeep = 0;
  int sumb = 0;
  uint64_t * ptbul = NULL;
  int i = 0; int j = 0;
  for (int h = 0; h <nvar; h++ )
  { ptbul = Bitdata[h];
    sumb = 0;
    for (i = 0; i < maxul;i++)
    {
      sumb += __builtin_popcountl(ptbul[i]);
    }
    sumvec[h]=sumb;
    if (sumb > support ) tokeep++;
  }
  Bdata = (uint64_t**) malloc(nvar*sizeof(uint64_t*));
  if (ordre == 'i' || ordre == 'd') 
  { int sumvecord [tokeep];
    for (i = 0; i < nvar;i++) {if (sumvec[i] > support) {Bdata[j] = Bitdata[i]; sumvecord[j] = sumvec[i]; varnames[j++] = varnames[i];}} 
    sortoutdata(Bdata, sumvecord,tokeep,varnames,ordre);
    return tokeep;
  }
  for (i = 0; i < nvar; i++) {if (sumvec[i]>support) {Bdata[j] = Bitdata[i]; varnames[j++] = varnames[i];} }
  
  return tokeep;
}

void creabitfield ( int *ind,int* mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbca)
{   auto * yb = &Ultab[0]; uint64_t u; int t = 0; auto * ybi = &Ultab[0];
auto * yd = &Ultab[64];
for (auto a =0; a < nbca ; a++)
{  
  for (yb = ybi; yb != yd; yb++)
  { 
    u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
  }
  nlgtabl[a] = u;
}
}

void creabitfieldr ( int *ind,int * mod,uint64_t * ptul, uint64_t * nlgtabl,  int nbdases,  int &reste)
{   
  auto * yb = &Ultab[0]; uint64_t u; int t = 0;
  auto * yd = &Ultab[64]; auto * ybi = &Ultab[0];
  for (auto a =0; a < nbdases ; a++)
  {  
    for (yb = ybi; yb != yd; yb++)
    { 
      u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;
    }
    nlgtabl[a] = u;
  }
  u = 0;
  for ( yb = ybi; yb != &Ultab[reste]; yb++) 
  {
    u ^= (-((ptul[ind[t]]>>mod[t++]) & 1UL) ^ u) & *yb;}
  nlgtabl[nbdases] = u;
  
}



void init_Rnodes (rnodes* outrnodes, int* vecsom,  int* tindices ,int * mod,uint64_t** Bitdata, int &nbleft, int& nbcases, int Mins,  int &reste)
{ 
  int sit = 0;
  int sum = 0; 
  if (reste ==0)
  {for (auto i=0; i < nbleft;i++)
  { sum = vecsom[i];
    if (sum<minsup) {outrnodes[i]= rnodes(0,0,NULL); sit++; continue;}
    nbfreq++;
    uint64_t * ptu = (uint64_t*) malloc(nbcases*sizeof(uint64_t));
    creabitfield(tindices,mod,Bitdata[sit],ptu,nbcases);
    outrnodes[i]= rnodes(sum,1,ptu);
    sit++;
  }
  }
  
  else {
    for (auto i=0; i < nbleft;i++)
    { sum = vecsom[i];
      if (sum<minsup) {outrnodes[i]= rnodes(0,0,NULL); sit++; continue;}
      nbfreq++;
      uint64_t * ptu = (uint64_t*) malloc(nbcases*sizeof(uint64_t));
      creabitfieldr(tindices,mod,Bitdata[sit],ptu,nbcases-1,reste);
      outrnodes[i]= rnodes(sum,1,ptu);
      sit++;
    }
  }
}

void depthwalk (pnodes * Tree, uint64_t * cand, rnodes* tabrnodes)
{ int elem = Tree->litem;
  bool d = !Tree->brother;
  if (!tabrnodes[elem].indic) {if (!d) {depthwalk(Tree->brother, cand, tabrnodes);} return;}
  int supcand = 0;
  auto * tib = tabrnodes[elem].tab;
  auto * ed = &cand[nbc];
  for (auto * tob = cand; tob != ed; ++tob,++tib )
  {supcand += __builtin_popcountl(*tob&*tib);}
  if (supcand<minsup) 
  {tabrnodes[elem].indic=0;
    if (!d) depthwalk(Tree->brother, cand, tabrnodes);
    tabrnodes[elem].indic=1;
    return;
  }          
  nbfreq++;
  if (!Tree->son)
  {Tree->son = new pnodes(supcand, curlitem,NULL);
    if (!d) {depthwalk(Tree->brother, cand, tabrnodes);}
    return;
  }   
  pnodes * trpnodes = Tree->son;
  uint64_t  * ptab = (uint64_t*) malloc (nbc*sizeof(uint64_t)) ;
  auto * pptab = &ptab[0];
  tib = tabrnodes[elem].tab;
  for (auto * tob = cand; tob != ed; ++tob,++tib,++pptab ){*pptab = *tob&*tib;}        
  depthwalk(trpnodes,ptab,tabrnodes);
  free(ptab);
  Tree->son = new pnodes (supcand,curlitem,trpnodes);
  if (!d) depthwalk(Tree->brother, cand, tabrnodes);
  return;       
}


void rootwalk (pnodes & root ,rnodes * rootnodes,  int nbcases, int nbleft )
{ pnodes * temp = root.brother;
  pnodes * trpn = NULL;
  for (int i =nbleft-1; i > -1 ; i--)
  {if (rootnodes[i].indic)
  {   trpn = temp->son;
    if ( trpn == NULL) temp->son = new pnodes (rootnodes[i].sup,curlitem,NULL);
    else 
    {  
      depthwalk(trpn,rootnodes[i].tab,rootnodes);
      temp->son = new pnodes (rootnodes[i].sup,curlitem,trpn);
    }
    
  }
  temp = temp->brother;
  }
}

void extract_and_erase_freq (pnodes*& Tree)
{ 
  if (Tree->son) {extract_and_erase_freq(Tree->son);}
  if (Tree->brother) {extract_and_erase_freq(Tree->brother);}
  delete Tree;
}


void extract_and_erase_freq (pnodes*& Tree, std::string str, std::string *  listenom, std::string * namevalue, int * supvalue, int* sizevalue, int* litemvalue , int size ,int& situeur )
{ 
 std::string strtemp = str + listenom[Tree->litem];
 namevalue[situeur] = strtemp;
 litemvalue[situeur] = Tree->litem;
 sizevalue [situeur] = size;
 supvalue[situeur++] = Tree->sup;
 
 
 if (Tree->son) {extract_and_erase_freq(Tree->son,strtemp,listenom,namevalue,supvalue,sizevalue,litemvalue,size+1,situeur);}
 if (Tree->brother) {extract_and_erase_freq(Tree->brother,str,listenom,namevalue,supvalue,sizevalue,litemvalue,size,situeur);}
delete Tree;

}


pnodes Bitprefrec (uint64_t** Bitdata, std::vector<std::string> &varnames, int maxul, int nvar, char ordre)
{  curlitem = 0; nbfreq = 0;
Rcpp::Rcout << "Support value is " << minsup << ". Checking unfrequent items ... \n ";
int sum = 0;
auto a = 0; int b= 0;
auto e = 0;
auto n = 0;
uint64_t ** Bdata = NULL ;
auto * pul = &Bitdata[0][0];
auto * pulb = &Bitdata[0][0];
auto * pulc = &Bitdata[0][0];
int keep = init_prefixtree(Bitdata, Bdata ,varnames, minsup,nvar,maxul,ordre);
Rcpp::Rcout <<"Done The Supportvalue is " << minsup << ". Start PrefRec with  " << keep << " frequents variables ..." <<  std::endl;
if(keep ==0)
{
  std::vector<std::string> outempty;
  outempty.push_back("empty");
  Rcpp::Rcout << "Not a single frequent itemSet, please try with a lower relative minSup value. " << std::endl;
  return pnodes (0,0,NULL);
}

clock_t stimer, etimer;
stimer = clock();
pnodes rootalpha (0,0,NULL);
for (int j=0; j<keep;j++)
{ 
  sum = 0;
  uint64_t u;
  curlitem = j;
  nbfreq++;
  pul = Bdata[j];
  pulc = &Bdata[j][maxul];
  for (a = 0; a < maxul;a++)
  {sum+= __builtin_popcountl(pul[a]);}
  nbc = sum/64;
  int reste = sum-64*nbc;
  if (reste) nbc ++;
  e = 0;
  rnodes tabrn [j+1]; 
  int indices [sum];
  int mod [sum];
  for (a=0; a < maxul; a++)
  { u = pul[a];
    for (n = 0; n < 64; n++)
      { if ((u >> n) & 1UL) {indices[e]= a; mod[e++]=n;} }
  }
  
  int rootsum [j+1]{0};
  int s = 0;
  for (a = 0; a < j; a++)
  {
    s = 0; pulb = Bdata[a];
    for (pul = Bdata[j]; pul != pulc; ++pul,++pulb)
    { s += __builtin_popcountl(*pul & *pulb);}
    rootsum[a] = s;
  }
  
  init_Rnodes(tabrn,rootsum,indices,mod,Bdata,j,nbc,minsup,reste);
  rootwalk(rootalpha,tabrn,nbc,j);
  pnodes * nouvnod = new pnodes (sum,curlitem,NULL);
  nouvnod->brother = rootalpha.brother;
  rootalpha.brother = nouvnod;
  
  for (a = 0; a < j; a++)
  {if (tabrn[a].tab) free(tabrn[a].tab);}
}

etimer = clock();
double diffe = double(etimer - stimer)/double(CLOCKS_PER_SEC);
Rcpp::Rcout << "End PrefRec, there is  " << nbfreq <<" sets, extracted int  " <<diffe << " secs"<< std::endl;
free (Bdata);

return rootalpha; 
}







// [[Rcpp::export]]
Rcpp::List Prefrec (std::vector<std::string> Data, double relativeSup, char deli, char ordre)
{
  
  int nrows = Data.size();
  int maxul = 0;
  int nvar = 0;
  uint64_t ** Bitdata = NULL;
  std::vector<std::string> varnames;
  for (int i = 0; i < 64; i++) {Ultab[i] = (1UL <<i);}
  
  Init_data(Data, nrows, maxul, nvar, Bitdata, varnames, deli );
  minsup = nrows*relativeSup;
  pnodes outnodes = Bitprefrec(Bitdata,varnames,maxul,nvar,ordre);
  if (nbfreq == 0)
    {
      std::vector<int> outnul(1);
      return Rcpp::List::create( Rcpp::Named("Abort")= outnul);
    }
  

 int choice = 1;
 if (nbfreq > 100000000)
   { Rcpp::Function readline = base["readline"];  
     Rcpp::Rcout << "The number of frequent is really large, be sure you have enought available memory before continue \n If you wan't to continue , press 1 else press 0 \n";
     choice = std::stoi(Rcpp::as<std::string>(readline(" Continue 1 , Stop 0 "))); Rcpp::Rcout << std::endl; }
 if (choice ==0 ) {extract_and_erase_freq(outnodes.brother); std::vector<int> outnul(1);  
     return Rcpp::List::create( Rcpp::Named("Abort")= outnul);}

 
 std::vector<int> supportvalue (nbfreq);
 std::vector<std::string> namevalue (nbfreq);
 std::vector<int> sizevalue (nbfreq);
 std::vector<int> litemvalue (nbfreq);
 int situ = 0;
 std::string startstring ;
 startstring.reserve (50);
 Rcpp::Rcout << "Clean and build R object" << std::endl;
 extract_and_erase_freq(outnodes.brother,startstring,&varnames[0],&namevalue[0],&supportvalue[0],&sizevalue[0],&litemvalue[0],1,situ);
 std::vector<double> relativesupportvalue (nbfreq);
 for (int q =0; q < nbfreq; q++)
 {relativesupportvalue[q] = (double)supportvalue[q]/(double)nrows;}
 
 DataFrame freqset = Rcpp::DataFrame::create(Rcpp::Named("id") = namevalue,
                                             Rcpp::Named("support") = supportvalue,
                                             Rcpp::Named("relative_support")=relativesupportvalue);
 
 DataFrame ruleindic = Rcpp::DataFrame::create(
                                              Rcpp::Named("Rsupport") = relativesupportvalue,
                                              Rcpp::Named("size") = sizevalue ,
                                              Rcpp::Named("litem") = litemvalue );
 
 
 return Rcpp::List::create( Rcpp::Named("frequent_itemset")= freqset,
                            Rcpp::Named("freqindic")=ruleindic,
                            Rcpp::Named("nbfreq")=nbfreq,
                            Rcpp::Named("item") = varnames );
 
}
