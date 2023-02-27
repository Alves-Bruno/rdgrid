#include<iostream>
#include<fstream>
#include<vector>
#include<sstream> // for ostringstream
#include<string>
#include<algorithm>
#include <cmath>
#include<cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <sys/stat.h>
#include <cstring>

#ifndef NORCPP
#include <Rcpp.h>
using namespace Rcpp;
#endif

unsigned int limit = 1;
using namespace std;

class point{
public:
  double x, y;

  template<class Archive>
  void serialize(Archive & archive)
  {
    archive(x, y); 
  }

  point(){
    x=0.0; y=0.0;
  }

  inline bool operator==(const point& p)
  {
    return ((this->x == p.x) && (this->y == p.y));
  }
};

class quad {
public: 
  quad **children;
  point a, b;
  vector<struct point> points;
  unsigned long N;
  bool has_children;

  template<class Archive>
  void serialize(Archive & archive)
  {
    archive(a, b, CEREAL_NVP(points), has_children, N); 
  }

  quad(){
    children = nullptr;
    N = -1;
    has_children = false;
  }

};

struct grid{
  quad *q0;
  vector<point> txrx;

  template<class Archive>
  void serialize(Archive & archive)
  {
    archive(CEREAL_NVP(txrx)); 
  }

  grid(){
    q0 = nullptr;
  }
};

std::string point_str(point *p){
  std::ostringstream pstr;
  // pstr << "(" << p->x << ", " << p->y << ")"; 
  pstr << p->x << ", " << p->y; 
  return(pstr.str());
}

void point_create(point *p, double x, double y){
  p->x=x; p->y=y;
}

void get_y(vector<double> &v, vector<point> &points){
  for(auto i : points){
    v.push_back(i.y);
  }
}
void get_x(vector<double> &v, vector<point> &points){
  for(auto i : points){
    v.push_back(i.x);
  }
}

void quad_print(quad *q){
  //cout << "Quad: " << endl;
  //cout << "\t" << point_str(&q->a) << " -> " << point_str(&q->b) << endl;
  cout << point_str(&q->a) << ", " << point_str(&q->b);
  cout << ", " << q->N << ", " << q->points.size() << ",  " << q->has_children << endl;
   //cout << q->npoints << endl;
  // if(q->points != nullptr){
  //   for(auto i : *(q->points)){
  //     //cout << i << endl;
  //     cout << "\t " << point_str((point*) &i) << endl;
  //   }
  // }
}

void print_all_sons(quad *q){
  quad_print(q);
  if(q->children!=nullptr){
    for(int i=0; i<4; i++)
      print_all_sons(q->children[i]);
  }
}


void quad_delete_all(quad *q){
  if(q->has_children){
    quad_delete_all(q->children[0]);
    quad_delete_all(q->children[1]);
    quad_delete_all(q->children[2]);
    quad_delete_all(q->children[3]);
    delete[] q->children;
  }

  //cout << "size: " << sizeof(q) << endl;
  delete(q);
  //delete(q->children);  
}

void grid_delete(grid *g){
  quad_delete_all(g->q0);
}

void quad_create(quad *q, point a, point b){
  q->a = a;
  q->b = b;
  q->children=nullptr;
  q->points=vector<point>();
  q->N=-1;
  q->has_children=false;
}

struct sort_by_x
{
    inline bool operator() (const point& p1, const point& p2)
    {
        return (p1.x < p2.x);
    }
};

struct sort_by_y
{
    inline bool operator() (const point& p1, const point& p2)
    {
        return (p1.y < p2.y);
    }
};


void split_sides(vector<point> &half_p, double mid, vector<point> &l, vector<point> &r){
  vector<double> half;
  get_x(half, half_p);
  auto xbound_top = lower_bound(half.begin(), half.end(), mid,
            [](const double& test, double value)
            {
                return test <= value;
            });
  int xbound_int = xbound_top - half.begin();
  int points_size = half.size();

  l.clear(); r.clear();
  l = vector<point>{&half_p[0], &half_p[xbound_int]};
  sort(l.begin(), l.end(), sort_by_y());
  r = vector<point>{&half_p[xbound_int], &half_p[points_size]};
  sort(r.begin(), r.end(), sort_by_y());
}

// Split the given q in 4 - once
void split_quad_once(grid *g, quad *q){

  // Create the 4 quads
  q->has_children=true;
  q->children = new(quad*[4]);
  for(int i=0; i<4;i++)
    q->children[i] = new(quad);
  
  point m, m_top, m_bot, m_left, m_right;
  point_create(&m, q->a.x + (q->b.x - q->a.x)/2.00, q->a.y + (q->b.y - q->a.y)/2.00);
  point_create(&m_top, m.x, q->b.y);
  point_create(&m_bot, m.x, q->a.y);
  point_create(&m_left, q->a.x, m.y);
  point_create(&m_right, q->b.x, m.y);

  quad_create(q->children[2], m_left, m_top);
  quad_create(q->children[3], m, q->b);
  quad_create(q->children[0], q->a, m);
  quad_create(q->children[1], m_bot, m_right);
  
  vector<double> y_values;
  vector<point> copy;
  for(auto i : q->points)
    copy.push_back(i);
  
  get_y(y_values, q->points);
  //for(auto i : y_values)
  //  cout << i << endl;
  auto ybound = lower_bound(y_values.begin(), y_values.end(), m.y,
            [](const double& test, double value)
            {
                return test <= value;
            });
  int ybound_int = ybound - y_values.begin();
  int points_size = q->points.size();

  vector<point> bot_half = {&(q->points)[0], &(q->points)[ybound_int]};
  sort(bot_half.begin(), bot_half.end(), sort_by_x());
  vector<point> top_half = {&(q->points)[ybound_int], &(q->points)[points_size]};
  sort(top_half.begin(), top_half.end(), sort_by_x());

  //  for(auto i:top_half)
  //  cout << "\t" << point_str(&i) << endl;
  
  //vector<point> top_half_l, top_half_r, bot_half_l, bot_half_r;
  split_sides(top_half, m.x, q->children[2]->points, q->children[3]->points);
  split_sides(bot_half, m.x, q->children[0]->points, q->children[1]->points);

  // Update N:
  for(int i=0; i<4; i++){
    q->children[i]->N = q->children[i]->points.size();
  }  
}

// Split the given q in 4 
void split_quad(grid *g, quad *q, bool keep_going){

  if(!keep_going){
    q->N = q->points.size();
    //g->result.push_back(q);
    //quad_print(q);
    return;
  }

  // Create the 4 quads
  q->has_children=true;
  q->children = new(quad*[4]);
  for(int i=0; i<4;i++)
    q->children[i] = new(quad);
  
  point m, m_top, m_bot, m_left, m_right;
  point_create(&m, q->a.x + (q->b.x - q->a.x)/2.00, q->a.y + (q->b.y - q->a.y)/2.00);
  point_create(&m_top, m.x, q->b.y);
  point_create(&m_bot, m.x, q->a.y);
  point_create(&m_left, q->a.x, m.y);
  point_create(&m_right, q->b.x, m.y);

  quad_create(q->children[2], m_left, m_top);
  quad_create(q->children[3], m, q->b);
  quad_create(q->children[0], q->a, m);
  quad_create(q->children[1], m_bot, m_right);
  
  vector<double> y_values;
  vector<point> copy;
  for(auto i : q->points)
    copy.push_back(i);
  
  get_y(y_values, q->points);
  //for(auto i : y_values)
  //  cout << i << endl;
  auto ybound = lower_bound(y_values.begin(), y_values.end(), m.y,
            [](const double& test, double value)
            {
                return test <= value;
            });
  int ybound_int = ybound - y_values.begin();
  int points_size = q->points.size();

  vector<point> bot_half = {&(q->points)[0], &(q->points)[ybound_int]};
  sort(bot_half.begin(), bot_half.end(), sort_by_x());
  vector<point> top_half = {&(q->points)[ybound_int], &(q->points)[points_size]};
  sort(top_half.begin(), top_half.end(), sort_by_x());

  //  for(auto i:top_half)
  //  cout << "\t" << point_str(&i) << endl;
  
  //vector<point> top_half_l, top_half_r, bot_half_l, bot_half_r;
  split_sides(top_half, m.x, q->children[2]->points, q->children[3]->points);
  split_sides(bot_half, m.x, q->children[0]->points, q->children[1]->points);

  // Update N:
  for(int i=0; i<4; i++){
    q->children[i]->N = q->children[i]->points.size();
  }

  if(keep_going){
    if(q->children[0]->points.size() > limit)
      split_quad(g, q->children[0], true);
    else split_quad(g, q->children[0], false);
           
    if(q->children[1]->points.size() > limit)
      split_quad(g, q->children[1], true);
    else split_quad(g, q->children[1], false);
  
    if(q->children[2]->points.size() > limit)
      split_quad(g, q->children[2], true);
    else split_quad(g, q->children[2], false);
  
    if(q->children[3]->points.size() > limit)
      split_quad(g, q->children[3], true);
    else split_quad(g, q->children[3], false);

  } 
  
} 

void grid_refine(grid *g){

  if(g->q0->points.size() > limit)
    split_quad(g, g->q0, true);

  //cout << "=================" << endl;
  //print_csv(g);
}

void walk_grid(quad *q, vector<quad*> &order){

  order.push_back(q);
  if(q->children != nullptr){
    q->has_children = true;
    walk_grid(q->children[0], order);
    walk_grid(q->children[1], order);
    walk_grid(q->children[2], order);
    walk_grid(q->children[3], order);
  } else {
    q->has_children = false;
  }
  
}

struct quad_queue {
  quad *q;
  int nsons;
};

void mount_quad_tree(grid *g, vector<quad*> quad_list){
  vector<quad_queue*> queue_fathers = {};
  
  for(auto i : quad_list){

    if(i==nullptr){
      cout << "MAYBE A BAD ALLOC - RUN!" << endl;
      exit(1);
    }
    // for(auto j : queue_fathers){
    //   cout << "queue: "; quad_print(j->q);
    // }

    if(i->has_children){
      if(queue_fathers.size() == 0){
        quad_queue *q = new(quad_queue);
        q->q=i; q->nsons=0;
        queue_fathers.push_back(q);
      }
      else {
        quad_queue *q = queue_fathers.back();
        //q.nsons+=1;

        if(q->q->children == nullptr)
          q->q->children = new(quad*[4]);

        q->q->children[q->nsons] = i;
        q->nsons = q->nsons + 1;

        if(q->nsons == 4){
          queue_fathers.pop_back();
          delete(q);
        }

        quad_queue *q2 = new(quad_queue);
        q2->nsons=0; q2->q=i;
        if(i->children == nullptr)
          i->children = new(quad*[4]);
        queue_fathers.push_back(q2);
        
      }
    } else {

      quad_queue *q = queue_fathers.back();

      if(q->q->children == nullptr)
        q->q->children = new(quad*[4]);

      q->q->children[q->nsons] = i;
      q->nsons = q->nsons + 1;

      if(q->nsons == 4){
        queue_fathers.pop_back();
        delete(q);
      }
        
    }
  }

  g->q0 = quad_list[0];
}

bool grid_import_struct(grid *g, const char *fname){

  std::ifstream is(fname);
  cereal::BinaryInputArchive ar(is);

  vector<quad*> quad_list;
  
  bool next = true;

  quad *q;
  while(next){
    try{
      try{
        q = new(quad);
      } catch (std::bad_alloc&) {
        // Handle error
        cout << "BAD ALLOC -------- RUUUUUNNN" << endl;
      }
      //q->points = vector<point>{};
      //quad_cereal qc;
      ar(*q);
      //quad_cereal_copy(q, qc);
      //quad_print(q);
      q->children = nullptr;
      quad_list.push_back(q);
    } catch(cereal::Exception cerror){
      //cout << "Stop!" << endl;
      next=false;
      delete(q);
      is.close();
    }
  }

  //for(auto i : quad_list){
  //  quad_print(i);
  //}
  
  if(quad_list.size() > 0) {
    mount_quad_tree(g, quad_list);
    //cout << "-----------------" << endl;
    //print_all_sons(g->q0);
    
    return true;
  } else return(false);
}

// Copy the points from the quad-tree leaves
// to the ones at the above sub-trees. You
// should consider leave_points_at_leaves,
// if you want to save memory ad disk space.
void rebuild_points(quad *q){
  if(q->children != nullptr){
    q->points.clear();
    for(int i=0; i<4; i++){
      rebuild_points(q->children[i]);
      q->points.insert(q->points.end(), q->children[i]->points.begin(), q->children[i]->points.end());
    }
    //remove(q->points);
  }
  q->N = q->points.size();
}

void leave_points_at_leaves(quad *q){
  if(q->children != nullptr){
    q->points.clear();
    q->N = 0;
    for(int i=0; i<4; i++){
      leave_points_at_leaves(q->children[i]);
      q->N = q->N + q->children[i]->N;
    }
  } else {
    q->N = q->points.size();
  }
}

void grid_export_struct(grid *g, const char *fname){
  // Write the data to the archive
  std::ofstream file(fname);
  cereal::BinaryOutputArchive archive(file);

  // Copy points from the leaves to the upper quads
  //rebuild_points(g->q0);

  // Leave points at the the leaves, just copy the
  // number of points to the upper quads.
  leave_points_at_leaves(g->q0);

  vector<quad*> grid_quads;
  walk_grid(g->q0, grid_quads);

  for(auto i : grid_quads)
    archive(*i);
      // quad_print(i);
  //  archive(*i);
  file.close();
  
}

void grid_export_txrx_points(grid *g, const char *fname){
  std::ofstream file(fname);
  cereal::BinaryOutputArchive archive(file);
  archive(*g);
  file.close();
}

void grid_import_txrx_points(grid *g, const char *fname){
  std::ifstream is(fname);
  cereal::BinaryInputArchive ar(is);
  ar(*g);
  is.close();
}

void create_common_grid(grid *A, grid *B, quad *qa, quad *qb){

  //quad_print(qa);
  //quad_print(qb);
  //cout << "--------" << endl;
  
  if(qa->has_children && qb->has_children){
    for(int i=0; i < 4; i++){
      create_common_grid(A, B,
                       qa->children[i],
                       qb->children[i]);
    }
  }

  else if(qa->has_children){
    split_quad_once(B, qb);
    create_common_grid(A, B,
                     qa,
                     qb);
      
  }

  else if(qb->has_children){
    split_quad_once(A, qa);
    create_common_grid(A, B,
                     qa,
                     qb);

  }


      
}

void grid_get_leaves(quad *q, vector<quad*> &leaves){

  if(q!=nullptr){
    if(q->children==nullptr){
      leaves.push_back(q);
    } else {
      grid_get_leaves(q->children[0], leaves);
      grid_get_leaves(q->children[1], leaves);
      grid_get_leaves(q->children[2], leaves);
      grid_get_leaves(q->children[3], leaves);
    }

  } else{ 
    cout << "Quad not allocated" << endl;
    exit(1);
  }
}

// As pointed by MultiRRomero at:
// https://stackoverflow.com/questions/5254838/calculating-distance-between-a-point-and-a-rectangular-box-nearest-point
double quad_dist_to_point(quad *q, point *p){
  vector<double> dx = {q->a.x - p->x, 0, p->x - q->b.x};
  vector<double> dy = {q->a.y - p->y, 0, p->y - q->b.y};

  double max_dx = *max_element(begin(dx), end(dx)); 
  double max_dy = *max_element(begin(dy), end(dy));
  
  return(sqrt(max_dx * max_dx + max_dy * max_dy));

}

void find_composed_N(grid *A, grid *B, quad *qa, quad *qb, quad *qc){
  //quad_print(qa);
  //quad_print(qb);
  //cout << endl;

  vector<double> A_dist, B_dist;
  //  quad_dist_to_point(A->txrx, )
  for(auto point_a : A->txrx){
    A_dist.push_back(quad_dist_to_point(qa, &point_a));
  }

  for(auto point_b : B->txrx){
    B_dist.push_back(quad_dist_to_point(qb, &point_b));
  }

  double min_a = *min_element(begin(A_dist), end(A_dist));
  double min_b = *min_element(begin(B_dist), end(B_dist));

  if(qa->points.size() >= qb->points.size()){
    qc->points = vector<point>{qa->points};
    qc->N = qa->points.size();
  }
  // B is more refined 
  else{
    qc->points = vector<point>{qb->points};
    qc->N = qb->points.size();
  }

  // if(min_a < min_b){
  //   qc->points = vector<point>{qa->points};
  //   qc->N = qa->points.size();
  // }
  // else if(min_a > min_b) {
  //   qc->points = vector<point>{qb->points};
  //   qc->N = qb->points.size();
  // }
  // else if(min_a == min_b){
  //   // A is more refined 
  //   if(qa->points.size() >= qb->points.size()){
  //     qc->points = vector<point>{qa->points};
  //     qc->N = qa->points.size();
  //   }
  //   // B is more refined 
  //   else{
  //     qc->points = vector<point>{qb->points};
  //     qc->N = qb->points.size();
  //   }
  // }
}

void copy_quad(quad *src, quad* dest){
  dest->has_children = src->has_children;
  dest->points = vector<point> {src->points};
  dest->N = src->points.size();
  dest->a = src->a;
  dest->b = src->b;
  if(src->children != nullptr){
    dest->children = new(quad*[4]);
    for(int i = 0; i < 4; i++){
      dest->children[i] = new(quad);
      copy_quad(src->children[i], dest->children[i]);
    }
    
  }
}

void copy_grid(grid *src, grid *dest){
  quad *d0 = new(quad);
  dest->q0 = d0;
  copy_quad(src->q0, dest->q0);
}

void remove(std::vector<point> &v)
{
    auto end = v.end();
    for (auto it = v.begin(); it != end; ++it) {
      end = std::remove(it + 1, end, *it);
    }
 
    v.erase(end, v.end());
}

void grid_compose(grid *A, grid *B, grid *C){

  quad *qa = A->q0;
  quad *qb = B->q0;

  vector<quad*> A_leaves, B_leaves, C_leaves;

  create_common_grid(A, B, qa, qb);
  copy_grid(A, C);
  quad *qc = C->q0;
  
  grid_get_leaves(qa, A_leaves);
  grid_get_leaves(qb, B_leaves);
  grid_get_leaves(qc, C_leaves);

  for(auto i=0; i < C_leaves.size(); i++){
    find_composed_N(A, B, A_leaves[i], B_leaves[i], C_leaves[i]);
  }

  // Copy txrx points from A and B to C
  for(auto i : A->txrx){
    C->txrx.push_back(i);
  }
  for(auto i : B->txrx){
    C->txrx.push_back(i);
  }
  // Remove duplicates
  remove(C->txrx);

}

// int main() {

//   // quad *q = new(quad);
//   // point a;
//   // point_create(&a, 0.0, 0.0);
//   // q->points = vector<point>{};
//   // q->points.push_back(a);
//   // delete(q);
  
//   grid *a = new(grid);
//   grid *b = new(grid);
//   grid *c = new(grid);

//   point txrx[4];
//   point_create(&txrx[0], 4.0, 4.5);
//   point_create(&txrx[1], 6.0, 3.0);
//   point_create(&txrx[2], 4.0, 4.5);
//   point_create(&txrx[3], 8.0, 3.0);

//   a->txrx.push_back(txrx[0]);
//   a->txrx.push_back(txrx[1]);
//   b->txrx.push_back(txrx[2]);
//   b->txrx.push_back(txrx[3]);
  
//   grid_import_struct(a, "a.dgrid");
//   //grid_import_txrx_points(a, "a.txrx");
  
//   grid_import_struct(b, "b.dgrid");
//   //  grid_import_txrx_points(b, "b.txrx");

//   grid_compose(b, a, c);

//   print_all_sons(c->q0);
//   vector<quad*> leaves;
//   grid_get_leaves(c->q0, leaves);
//   for(auto i : leaves){
//     for(auto j : i->points)
//       cout << point_str(&j) << endl;
//   }
  
//   // print_all_sons(b->q0);
//   // for(auto i : b->q0->points){
//   //    cout << point_str(&i) << endl;
//   // }
  
//   grid_delete(a);
//   delete(a);
//   grid_delete(b);
//   delete(b);
//   grid_delete(c);
//   delete(c);
  
  
//   return 0;
// }

#ifndef NORCPP
// [[Rcpp::export]]
void dgrid_create(unsigned long int max_N,
                  double xmin, double ymin, double xmax, double ymax,
                  NumericVector txrx_x, NumericVector txrx_y,
                  NumericVector midpoints_x, NumericVector midpoints_y,
                  String fname){

  limit = max_N;
  point a, b;
  point_create(&a, (double) xmin, (double) ymin);
  point_create(&b, (double) xmax, (double) ymax);
  
  quad *master = new(quad);
  quad_create(master, a, b);

  grid *g = new(grid);
  g->q0=master;

  for(int i=0; i < txrx_x.size(); i++){
    point p;
    point_create(&p, (double) txrx_x[i], (double) txrx_y[i]);
    g->txrx.push_back(p);
  }
  
  for(int i=0; i < midpoints_x.size(); i++){
    point p;
    point_create(&p, (double) midpoints_x[i], (double) midpoints_y[i]);
    master->points.push_back(p);
  }

  sort(master->points.begin(), master->points.end(), sort_by_y());
  grid_refine(g);

  leave_points_at_leaves(g->q0);
  //print_all_sons(g->q0);
  
  fname.replace_last("$", ".rdgrid");
  grid_export_struct(g, fname.get_cstring());
  fname.replace_last("$", ".txrx");
  grid_export_txrx_points(g, fname.get_cstring());
  //cout << fpath.get_cstring() << fname.get_cstring() << endl; 
  //print_csv(&g);

  grid_delete(g);
  delete(g);
}

// [[Rcpp::export]]
void dgrid_compose(String Afname, String Bfname, String composed){
  grid *A = new(grid);
  grid *B = new(grid);
  grid *C = new(grid);

  Afname.replace_last("$", ".rdgrid");
  grid_import_struct(A, Afname.get_cstring());
  Afname.replace_last("$", ".txrx");
  grid_import_txrx_points(A, Afname.get_cstring());

  Bfname.replace_last("$", ".rdgrid");
  grid_import_struct(B, Bfname.get_cstring());
  Bfname.replace_last("$", ".txrx");
  grid_import_txrx_points(B, Bfname.get_cstring());

  grid_compose(A, B, C);

  composed.replace_last("$", ".rdgrid");
  grid_export_struct(C, composed.get_cstring());
  composed.replace_last("$", ".txrx");
  grid_export_txrx_points(C, composed.get_cstring());
  
  grid_delete(A); grid_delete(B);
  grid_delete(C);
  delete(A); delete(B), delete(C);
}

// [[Rcpp::export]]
NumericMatrix dgrid2table(String fname){
  
  grid *g = new(grid);
  fname.replace_last("$", ".rdgrid");
  grid_import_struct(g, fname.get_cstring());
  fname.replace_last("$", ".txrx");
  grid_import_txrx_points(g, fname.get_cstring());

  //print_all_sons(g->q0);
  
  vector<quad*> leaves;
  grid_get_leaves(g->q0, leaves);

  vector<point> points;
  for(auto i : leaves){
    for(auto j : i->points)
      points.push_back(j);
  }

  unsigned long int nrows = leaves.size() + points.size() + g->txrx.size();
  unsigned long int ncols = 6;
  NumericMatrix df(nrows, ncols);
  // Copy leaves quads
  for(int i = 0; i < leaves.size(); i++){
    df[0 * nrows + i] = leaves[i]->a.x;
    df[1 * nrows + i] = leaves[i]->a.y;
    df[2 * nrows + i] = leaves[i]->b.x;
    df[3 * nrows + i] = leaves[i]->b.y;
    df[4 * nrows + i] = leaves[i]->N;
    df[5 * nrows + i] = 0.0;
  }

  // Copy leaves points
  for(int i = leaves.size(); i < leaves.size() + points.size(); i++){
    df[0 * nrows + i] = points[i - leaves.size()].x;
    df[1 * nrows + i] = points[i - leaves.size()].y;
    df[2 * nrows + i] = 0.0;
    df[3 * nrows + i] = 0.0;
    df[4 * nrows + i] = 0.0;
    df[5 * nrows + i] = 1.0;
  }

  // Copy rxtx points
  for(int i = leaves.size() + points.size(); i < nrows; i++){
    df[0 * nrows + i] = g->txrx[i - leaves.size() - points.size()].x;
    df[1 * nrows + i] = g->txrx[i - leaves.size() - points.size()].y;
    df[2 * nrows + i] = 0.0;
    df[3 * nrows + i] = 0.0;
    df[4 * nrows + i] = 0.0;
    df[5 * nrows + i] = 2.0;
  }

  colnames(df) = CharacterVector::create("xmin", "ymin", "xmax", "ymax", "N", "type");

  //print_all_sons(g->q0);
  grid_delete(g);
  delete(g);
  return(df);
 
}

// [[Rcpp::export]]
unsigned long int dgrid_depth(String fname){
  grid *g = new(grid);
  fname.replace_last("$", ".rdgrid");
  grid_import_struct(g, fname.get_cstring());
  unsigned long int depth = g->q0->N;
  grid_delete(g);
  delete(g);
  return(depth);
}
#endif

#ifdef NORCPP

void dgrid_compose(char* Afname, char* Bfname, char* composed){
  grid *A = new(grid);
  grid *B = new(grid);
  grid *C = new(grid);

  vector<string> fnames;
  fnames.push_back(string{Afname} + ".rdgrid");
  fnames.push_back(string{Afname} + ".rdgrid.txrx");
  fnames.push_back(string{Bfname} + ".rdgrid");
  fnames.push_back(string{Bfname} + ".rdgrid.txrx");
  fnames.push_back(string{composed} + ".rdgrid");
  fnames.push_back(string{composed} + ".rdgrid.txrx");

  grid_import_struct(A, fnames[0].c_str());
  grid_import_txrx_points(A, fnames[1].c_str());

  grid_import_struct(B, fnames[2].c_str());
  grid_import_txrx_points(B, fnames[3].c_str());

  grid_compose(A, B, C);

  grid_export_struct(C, fnames[4].c_str());
  grid_export_txrx_points(C, fnames[5].c_str());
  
  grid_delete(A); grid_delete(B);
  grid_delete(C);
  delete(A); delete(B), delete(C);
}

unsigned long int dgrid_depth(char* fname){
  string sfname = string{fname};
  sfname = sfname + ".rdgrid";

  //grid_import_struct(g, sfname.c_str());

  std::ifstream is(sfname.c_str());
  cereal::BinaryInputArchive ar(is);

  bool next = true;
  quad *q;
    try{
      try{
        q = new(quad);
      } catch (std::bad_alloc&) {
        // Handle error
        cout << "BAD ALLOC -------- RUUUUUNNN" << endl;
      }
      //q->points = vector<point>{};
      //quad_cereal qc;
      ar(*q);
      //quad_cereal_copy(q, qc);
      //quad_print(q);
      q->children = nullptr;
    } catch(cereal::Exception cerror){
      //cout << "Stop!" << endl;
      next=false;
      delete(q);
      is.close();
    }
  
  unsigned long int depth = q->N;
  delete(q);
  return(depth);
}

int main(int argc, char* argv[]){

  //dgrid_depth
  //dgrid_compose
  // Check the number of parameters
  if (argc < 3 || argc > 6) {
    std::cerr << "Usage: " << endl;
    std::cerr << "\tCompose: " << argv[0] << " -c gridA gridB out" << std::endl;
    std::cerr << "\tGet depth: " << argv[0] << " -d grid" << std::endl;
    return 1;
  }

  if (argc == 3) {
    if(strcmp(argv[1], "-d") == 0){
      struct stat buffer;
      vector<string> in_files;
      in_files.push_back(string(argv[2]) + ".rdgrid");
      for(auto f : in_files){
        if(!stat(f.c_str(), &buffer) == 0){
          cerr << "File not found: " << f << endl;
          return 1;
        }
      }
      cout << dgrid_depth(argv[2]) << endl;
      
      
    }
    else{
      cerr << "Option not known: " << argv[1] << endl;
    }
  }

  if (argc == 5) {
    if(strcmp(argv[1], "-c") == 0){
      struct stat buffer;
      vector<string> in_files;
      in_files.push_back(string(argv[2]) + ".rdgrid");
      in_files.push_back(string(argv[2]) + ".rdgrid.txrx");
      in_files.push_back(string(argv[3]) + ".rdgrid");
      in_files.push_back(string(argv[3]) + ".rdgrid.txrx");

      for(auto f : in_files){
        if(!stat(f.c_str(), &buffer) == 0){
          cerr << "File not found: " << f << endl;
          return 1;
        }
      }
      dgrid_compose(argv[2], argv[3], argv[4]);
      
    }
    else{
      cerr << "Option not known: " << argv[1] << endl;
    }
  }
  return 0;
}
#endif
