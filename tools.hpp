#pragma once

#include <iostream>

using namespace std;

class stopwatch
{
  struct timespec start_time, stop_time;
  double elapsed;
  bool running = false;
  bool onscreen = false;
public:
  void click()
  {
    if(running)
    {
      clock_gettime(CLOCK_MONOTONIC, &(this->stop_time));
      this->running = false;
      this->onscreen = true;
    }
    else
    {
      clock_gettime(CLOCK_MONOTONIC, &(this->start_time));
      this->running = true;
      this->onscreen = false;
    }
  }
  double check()
  {
    if(onscreen)
    {
      this->elapsed = (this->stop_time.tv_sec - this->start_time.tv_sec);
      this->elapsed += (this->stop_time.tv_nsec - this->start_time.tv_nsec) / 1000000000.0;
      return this->elapsed*1000;
    }
    else
      return 0.0;
  }
  void print_time()
  {
    std::cout<<this->check()<<"ms\n";
  }
};

void pvect(vector<float> v)
{
	for (auto entry:v)
	{
		cout<<entry<<" ";
	}
	cout<<endl;
}

void parray(float* v, int N)
{
	for (int i =0; i<N; i++)
	{
		cout<<v[i]<<" ";
	}
	cout<<endl;
}

void printmat(float* v, int rows, int cols)
{
	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < cols; col++)
		{
			cout<<v[col+row*cols]<<" ";
		}
		cout<<"\n";
	}
}

int rc2ii(int row, int col, int N)
{
    return row*N + col;
}
