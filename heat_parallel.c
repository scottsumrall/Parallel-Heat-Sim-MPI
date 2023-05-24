/*
 * 2D Heat Transfer - Discrete Event Simulation
 * 
 * Scott Sumrall
 * Reaches steady state on step 179 (180 steps) given 0.001
 *
 * Problem parallelized by dividing grid into rows based on number of threads. Picked
 * this method for simplicity. Program will distribute rows to processes and account for
 * to uneven distribution of rows per thread. Situations with an uneven distribution may
 * be less optimal since individual rows are not divided but should be effective at
 * reasonable scale. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "gfx.h"
#include <mpi.h>


#define SIZE 300

#define HEATER_COUNT 2
#define HEATER_SIZE 50
#define MAX_ITERATIONS 10000
#define COLD 0.0
#define HOT  1.0
#define AMBIENT 0.5

#define K 0.25

float grid[SIZE][SIZE];
float new_grid[SIZE][SIZE];

// heaters
int heater_location_x[HEATER_COUNT];
int heater_location_y[HEATER_COUNT];

void initialize_walls() {
    for (int i = 0; i < SIZE; i++) {
        grid[i][0] = COLD;       // left wall
        grid[0][i] = COLD;       // top wall
        grid[i][SIZE-1] = COLD;  // right wall
        grid[SIZE-1][i] = COLD;  // bottom wall

        new_grid[i][0] = COLD;       // left wall
        new_grid[0][i] = COLD;       // top wall
        new_grid[i][SIZE-1] = COLD;  // right wall
        new_grid[SIZE-1][i] = COLD;  // bottom wall
    }
}

void initialize_room() {
    for (int i = 1; i < SIZE - 1; i++) {
        for (int j = 1; j < SIZE - 1; j++) {
            grid[i][j] = AMBIENT;
        }
    }
}

void drop_heaters() {
    srand(time(NULL));
    
    for (int i = 0; i < HEATER_COUNT; i++) {
        heater_location_x[i] = (rand() % (SIZE - (HEATER_SIZE + 1))) + 1;
        heater_location_y[i] = (rand() % (SIZE - (HEATER_SIZE + 1))) + 1;

        // heater location (top-left corner)
        for (int x = 0; x < HEATER_SIZE; x++) {
            for (int y = 0; y < HEATER_SIZE; y++) {
                grid[heater_location_x[i] + x][heater_location_y[i] + y] = HOT;
            }
        }
    }
}

void update_heaters() {    
    for (int i = 0; i < HEATER_COUNT; i++) {
        // heater location (top-left corner)
        for (int x = 0; x < HEATER_SIZE; x++) {
            for (int y = 0; y < HEATER_SIZE; y++) {
                new_grid[heater_location_x[i] + x][heater_location_y[i] + y] = HOT;
            }
        }
    }
}

void print_room() {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            printf("%.1f", grid[i][j]);
        }
        printf("\n");
    }
}

void display_room() {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            float intensity = 255 * grid[i][j];

            float red = intensity;
            float green = 0;               
            float blue = 255 - intensity;  

            gfx_color(red, green, blue);
            gfx_point(i, j);
        }
    }
    gfx_flush();
}

// Debugging
void print_grid() {
    // Print column headers
    printf("   ");
    for (int j = 0; j < SIZE; j++) {
        printf("%-4d", j);
    }
    printf("\n");

    // Print grid rows
    for (int i = 0; i < SIZE; i++) {
        printf("%-4d", i);  // Print row header

        for (int j = 0; j < SIZE; j++) {
            printf("%.1f ", grid[i][j]);
        }
        printf("\n");
    }
}

void copy_room() {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            grid[i][j] = new_grid[i][j];
        }
    }
}

int main(int argc, char **argv) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    float threshold = atof(argv[1]);

    // Determine the size of the grid for each process
    // Take into account uneven distribution
    int chunk_size = SIZE / size;
    int remainder = SIZE % size;

    int start_row, end_row;

    if (rank < remainder) {
        start_row = rank * (chunk_size + 1);
        end_row = start_row + chunk_size + 1;
    } else {
        start_row = rank * chunk_size + remainder;
        end_row = start_row + chunk_size;
    }

    initialize_walls();
    initialize_room();
    
    drop_heaters();

    if(rank == 0)
    {
        gfx_open(SIZE, SIZE, "Heat Transfer");
    }

    int step = 0;
    int isStable = 0;

    while (!isStable) 
    { 
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == 0)
        {
            display_room();
            printf("STEP %d\n", step);
        }
		
		//update grid
        for (int i = start_row -1; i < end_row; i++) {
            for (int j = 1; j < SIZE - 1; j++) {
                new_grid[i][j] = grid[i][j] + K * (grid[i + 1][j] + grid[i - 1][j] + grid[i][j + 1] + grid[i][j - 1] - 4 * grid[i][j]);
            }
        }
		
		//distribute updated grid to processes
        MPI_Allgather(&new_grid[start_row][0], chunk_size * SIZE, MPI_FLOAT, &new_grid[0][0], chunk_size * SIZE, MPI_FLOAT, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        
        initialize_walls();
        update_heaters();

        MPI_Barrier(MPI_COMM_WORLD);

        // stability check
        isStable = 1;
        for (int i = 0; i < SIZE; i++) {
            for (int j = 0; j < SIZE; j++) {
                float diff = fabs(grid[i][j] - new_grid[i][j]);
                if (diff > threshold) {
                    isStable = 0;
                    break;
                }
            }
            if (!isStable) break;
        }

        copy_room();
        step++;
    }

    // Clean up MPI
    MPI_Finalize();
   
    if (rank == 0) {
        while (1) {
            char c = gfx_wait();
            if (c == 'q') break;
        }
    }
}






