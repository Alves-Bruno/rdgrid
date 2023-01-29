#include<iostream>
#include<fstream>
#include<vector>
#include<sstream> // for ostringstream
#include<string>
#include<algorithm>
#include<cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>


//#include <Rcpp.h>
//using namespace Rcpp;

unsigned int limit = 1;
using namespace std;

struct point{
  double x, y;

  template<class Archive>
  void serialize(Archive & archive)
  {
    archive(x, y); 
  }
};

struct quad {
  quad **children;
  point a, b;
  vector<struct point> points;
  unsigned long N;
  bool has_children;

  template<class Archive>
  void serialize(Archive & archive)
  {
    archive(a, b, CEREAL_NVP(points), has_children); 
  }

};

struct grid{
  quad *q0;
  vector<quad*> result;
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
  cout << ", " << q->points.size() << ",  " << q->has_children << endl;
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
  if(q->children != nullptr){
    quad_delete_all(q->children[0]);
    quad_delete_all(q->children[1]);
    quad_delete_all(q->children[2]);
    quad_delete_all(q->children[3]);
    free(q->children);
  }

  free(q);  
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

// Split the given q in 4 
void split_quad(grid *g, quad *q, bool keep_going){

  if(!keep_going){
    q->N = q->points.size();
    g->result.push_back(q);
    //quad_print(q);
    return;
  }

  // Create the 4 quads
  q->has_children=true;
  q->children = (quad **) malloc(sizeof(quad) * 4);
  for(int i=0; i<4;i++)
    q->children[i] = (quad *) malloc(sizeof(quad));
  
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

  //cout << "size: " << bot_half_l.size() << endl; 

  // Update structs
  //q->children[0]->points = &bot_half_l;
  //q->children[1]->points = &bot_half_r;
  //q->children[2]->points = &top_half_l;
  //q->children[3]->points = &top_half_r;

  
  //(*q->points).clear();
  //q->points=nullptr;

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

void print_csv(grid *g){

  for(auto i : g->result)
    quad_print(i);
}
void grid_refine(grid *g){

  if(g->q0->points.size() > limit)
    split_quad(g, g->q0, true);

  //cout << "=================" << endl;
  //print_csv(g);
}

// [[Rcpp::export]]
// NumericMatrix dynamic_grid(
//                            unsigned long int max_N,
//                            double xmin, double ymin, double xmax, double ymax,
//                            NumericVector midpoints_x, NumericVector midpoints_y
// ) {

//   limit = max_N;
//   point a, b;
//   point_create(&a, (double) xmin, (double) ymin);
//   point_create(&b, (double) xmax, (double) ymax);
  
//   quad master;
//   quad_create(&master, a, b);

//   grid g;
//   g.q0=&master;

//   vector<point> points;
//   for(int i=0; i < midpoints_x.size(); i++){
//     point p;
//     point_create(&p, (double) midpoints_x[i], (double) midpoints_y[i]);
//     points.push_back(p);
//   }

//   sort(points.begin(), points.end(), sort_by_y());
//   master.points=&points;
//   grid_refine(&g);
//   //print_csv(&g);

//   unsigned long int nrows = g.result.size();
//   unsigned long int ncols = 5;
    
//   NumericMatrix df(nrows, ncols);
//   for(int i = 0; i < nrows; i++){
//     df[0 * nrows + i] = g.result[i]->a.x;
//     df[1 * nrows + i] = g.result[i]->a.y;
//     df[2 * nrows + i] = g.result[i]->b.x;
//     df[3 * nrows + i] = g.result[i]->b.y;
//     df[4 * nrows + i] = g.result[i]->N;
//   }

//   //df[0] = 1000;
//   //df[1] = -1000;

//   return(df);
// }

void walk_grid(quad *q, vector<quad*> &order){

  order.push_back(q);
  if(q->children != nullptr){
    walk_grid(q->children[0], order);
    walk_grid(q->children[1], order);
    walk_grid(q->children[2], order);
    walk_grid(q->children[3], order);
  } 
  
}

struct quad_queue {
  quad * q;
  int nsons;
};

void mount_quad_tree(grid *g, vector<quad*> quad_list){
  vector<quad_queue*> queue_fathers = {};
  
  for(auto i : quad_list){

    // for(auto j : queue_fathers){
    //   cout << "queue: "; quad_print(j->q);
    // }

    if(i->has_children){
      if(queue_fathers.size() == 0){
        quad_queue *q = (quad_queue*) malloc(sizeof(quad_queue));
        q->q=i; q->nsons=0;
        queue_fathers.push_back(q);
      }
      else {
        quad_queue *q = queue_fathers.back();
        //q.nsons+=1;

        if(q->q->children == nullptr)
          q->q->children = (quad**) malloc(sizeof(quad)*4);

        q->q->children[q->nsons] = i;
        q->nsons = q->nsons + 1;

        quad_queue *q2 = (quad_queue*) malloc(sizeof(quad_queue));;
        q2->nsons=0; q2->q=i;
        if(i->children == nullptr)
          i->children = (quad**) malloc(sizeof(quad)*4);
        queue_fathers.push_back(q2);
        
      }
    } else {

      quad_queue *q = queue_fathers.back();
      
      q->q->children[q->nsons] = i;
      q->nsons = q->nsons + 1;

      if(q->nsons == 4){
        free(q);
        queue_fathers.pop_back();
      }
        
    }

    //   else{
    //     quad *father = queue_fathers.back();
    //     if(father->children==nullptr){
    //       father->children = (quad**) malloc(sizeof(quad)*4);
    //     }
    //     father->children[nsons.back()] = i;
    //     nsons[nsons.size()-1] = nsons.back() + 1;

    //     queue_fathers.push_back(i);
    //     nsons.push_back(0);
    //   }
    // } else{
     
    //   quad *father = queue_fathers.back();
    //   if(father->children==nullptr){
    //     father->children = (quad**) malloc(sizeof(quad)*4);
    //   }
    //   father->children[nsons.back()] = i;
    //   nsons[nsons.size()-1] = nsons.back() + 1;
    // }
      
  }

  g->q0 = quad_list[0];
}

bool grid_import_struct(grid *g){

  std::ifstream is("out.bin");
  cereal::BinaryInputArchive ar(is);

  vector<quad*> quad_list;
  
  bool next = true;
  
  while(next){
    try{
      quad *q = (quad*) malloc(sizeof(quad));
      ar(*q);
      quad_print(q);
      q->children = nullptr;
      quad_list.push_back(q);
    } catch(cereal::Exception cerror){
      //cout << "Stop!" << endl;
      next=false;      
    }
  }

  if(quad_list.size() > 0) {
    mount_quad_tree(g, quad_list);
    cout << "-----------------" << endl;
    print_all_sons(g->q0);
    
    return true;
  } else return(false);
}


void grid_export_struct(grid *g){
  // Write the data to the archive
  std::ofstream file( "out.bin" );
  cereal::BinaryOutputArchive archive(file);

  vector<quad*> grid_quads;
  walk_grid(g->q0, grid_quads);

  for(auto i : grid_quads)
    archive(*i);
      // quad_print(i);
  //  archive(*i);
  
  
}

int main() {
  point a, b;
  point_create(&a, 0.00, 0.00);
  point_create(&b, 10.00, 10.00);
  
  quad *master;
  master = (quad *) malloc(sizeof(quad));
  quad_create(master, a, b);

  grid g;
  g.q0=master;

  point p[6];
  point_create(&p[0], 1.0, 4.0);
  point_create(&p[1], 2.0, 1.0);
  point_create(&p[2], 4.0, 4.0);
  point_create(&p[3], 3.0, 7.0);
  point_create(&p[4], 7.0, 6.0);
  point_create(&p[5], 6.0, 9.0);

  //vector<point> points;
  master->points.push_back(p[0]);
  master->points.push_back(p[1]);
  master->points.push_back(p[2]);
  master->points.push_back(p[3]);
  master->points.push_back(p[4]);
  master->points.push_back(p[5]);

  sort(master->points.begin(), master->points.end(), sort_by_y());

  // quad_print(&master);
  //point_print(p[0]);

  grid_refine(&g);
  //quad_delete_all(master);
  grid_export_struct(&g);

  
  //  std::stringstream ss;
  //  cereal::BinaryOutputArchive oarchive(ss); 
  // archive(p[0]);
  
    
  grid_delete(&g);
  //print_csv(&g);

  grid newg;
  grid_import_struct(&newg);
  grid_delete(&newg);
  
  return 0;
}

