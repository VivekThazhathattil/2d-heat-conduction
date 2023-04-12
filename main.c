/* 2D Transient Heat Equation for steel plate solver via Finite-Difference scheme */

#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <time.h>

#define BLK "\e[0;30m"
#define RED "\e[0;31m"
#define GRN "\e[0;32m"
#define YEL "\e[0;33m"
#define BLU "\e[0;34m"
#define MAG "\e[0;35m"
#define CYN "\e[0;36m"
#define WHT "\e[0;37m"

#define gotoxy(x,y) printf("\033[%d;%dH", (y), (x))

typedef struct TEMP_S{
  int x, y;
  int num_snapshots;
  double* series;
} temp_t;

typedef struct WINSIZE_S{
  int rows, cols;
} winsize_t;


temp_t* initialize_temperature_series(int x, int y, int n){
  temp_t* temp = (temp_t*) malloc(sizeof(temp_t));
  temp->x = x;
  temp->y = y;
  temp->num_snapshots = n;
  temp->series = (double*) calloc(n, sizeof(double));
  return temp;
}

void delete_temperature_series(temp_t* t){
  free(t->series);
  free(t);
}

temp_t*** initialize_temperature_field(int nx, int ny, int n_snapshots){
  temp_t*** temp_field = (temp_t***) malloc(sizeof(temp_t**)*nx);

  for(int i = 0; i < nx; ++i){
    temp_field[i] = (temp_t**) malloc(sizeof(temp_t*)*ny);
    for(int j = 0; j < ny; ++j)
      temp_field[i][j] = initialize_temperature_series(i, j, n_snapshots);
  }
  return temp_field; 
}

void delete_temperature_field(int nx, int ny, temp_t*** t){
  for(int i = 0; i < nx; ++i){
    for(int j = 0; j < ny; ++j)
      delete_temperature_series(t[i][j]);
    free(t[i]);
  }
  free(t);
}

double* linspace(int lo, int hi, int n){
  if(n <= 1)
    return NULL;
  double* vals = (double*) malloc(sizeof(double)*n);
  double step = (double)(hi - lo)/ (n - 1);
  for(int i = 0; i < n; ++i)
    vals[i] = lo + step*i;
  return vals;
}

char* get_specific_color(double num){
  if(0 <= num && num <= 273)
    return BLK;
  else if(273 < num && num <= 300)
    return BLU;
  else if(300 < num && num <= 350)
    return CYN;
  else if(350 < num && num <= 400)
    return WHT;
  else if(400 < num && num <= 450)
    return GRN;
  else if(450 < num && num <= 500)
    return YEL;
  else if(500 < num && num <= 550)
    return RED;
  else if(550 < num && num <= 13000)
    return MAG;
}

winsize_t get_win_size(void){
  winsize_t ws;
  struct winsize win;     
  ioctl(0, TIOCGWINSZ, &win);
  ws.rows = win.ws_row;
  ws.cols = win.ws_col;
  return ws;
}

void draw_cells(char* c, temp_t*** temp_field, int nx, int ny, int ti){
  for(int i = 0; i < ny; ++i){
    for(int j = 0; j < nx; ++j){
      char* color = get_specific_color(temp_field[j][i]->series[ti]);
      gotoxy(j + 1, i + 1);
      printf("%s%s", color, c);
    }
  }
}

int main(){
  /* Physical params */
  double k = 1.172e-5; // for steel with 1% carbon
  double lx = 0.5; // length
  double ly = 0.5; // width

  /* Numerical params */
  winsize_t win = get_win_size();
  int nx = win.cols; // num points in x direction
  int ny = win.rows; // num points in y direction
  printf("%d, %d\n", nx, ny);
  double dt = 0.1; // time step
  double tf = 1000.0; // final time
                    //
  /* Draw param */
  //int frame_update_delay = 4096*8;
  int frame_update_delay = 0;
  char* blk_sym = "■";
  //char* blk_sym = "█";
           
  /* Boundary conditions (Dirichlet) */
  int temp0 = 273; // K at time T = 0
  int temp1 = 273; // Top Boundary
  int temp2 = 273; // Bottom Boundary
  int temp3 = 1000; // Left Boundary
  int temp4 = 273; // Right Boundary

  /* Cell dimensions */
  double dx = lx/nx;
  double dy = ly/ny;

  /* Courant numbers */
  double r1 = k*dt/(dx*dx);
  double r2 = k*dt/(dy*dy);

  if(r1 > 0.5 || r2 > 0.5){
    printf("Courant number error: Unstable solution. Exiting...");
    exit(1);
  }

  int num_tsteps = (int)(tf/dt);

  /* initialize temperature field */
  temp_t*** temp_field = initialize_temperature_field(nx, ny, num_tsteps);

  /* set initial condition */
  for(int i = 0; i < nx - 1; ++i)
    for(int j = 0; j < ny - 1; ++j)
      temp_field[i][j]->series[0] = temp0;

  /* set boundary conditions */

  /* Top and Bottom */
  for(int i = 0; i < nx; ++i){
    for(int t = 0; t < num_tsteps; ++t){
      temp_field[i][0]->series[t] = temp1;
      temp_field[i][ny-1]->series[t] = temp2;
    }
  }
  
  /* Left and Right */
  for(int j = 0; j < ny; ++j){
    for(int t = 0; t < num_tsteps; ++t){
      temp_field[0][j]->series[t] = temp3;
      temp_field[nx-1][j]->series[t] = temp4;
    }
  }

  /* Generate 2d mesh */
  //double* line_x = linspace(0, lx, nx);
  //double* line_y = linspace(0, ly, ny);
  //double** mesh_x = generate_mesh(line_x, nx, ny, 0);
  //double** mesh_y = generate_mesh(line_y, nx, ny, 1);
  
  system("clear");
  /* Main time-loop */
  for(int t = 0; t < num_tsteps; ++t){
    for(int i = 1; i < nx-1; ++i){
      for(int j = 1; j < ny-1; ++j){
        double d2_temp_dx2 = (temp_field[i+1][j]->series[t] 
          - 2*temp_field[i][j]->series[t] 
          + temp_field[i-1][j]->series[t])/(dx*dx);
        
        double d2_temp_dy2 = (temp_field[i][j+1]->series[t] 
          - 2*temp_field[i][j]->series[t] 
          + temp_field[i][j-1]->series[t])/(dy*dy);

        temp_field[i][j]->series[t+1] = k*dt*(d2_temp_dx2 + d2_temp_dy2) + temp_field[i][j]->series[t];
      }
    }
    draw_cells(blk_sym, temp_field, nx, ny, t);
    usleep(frame_update_delay);
  }

 // for(int t = 0; t < temp_field[1][2]->num_snapshots; ++t)
 //   printf("%f\n", temp_field[1][2]->series[t]);

  /* free temperature field and associated series */
  delete_temperature_field(nx, ny, temp_field);
  /* free other dynamically allocated variables */
  //free(line_x);
  //free(line_y);
  return 0;
}
