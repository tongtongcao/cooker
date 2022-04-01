//void findMaximum(double *array, int size, int *index, double *max);
void findMaximum(vector<double>array, int size, int *index, double *max);
void kera()
  
{
  //const int size=7;
  int index; double max;
  // double array[size]={1,2,9,7,3,15}; //= new double[size];
  // double *array = new double[7];
  vector<double>array;
  
  for(int j=0;j<7;j++)
    {
      // array[j]=500.0-j;
      array.push_back(500-j);
      cout<<array[j]<<endl;
    }
  // int size=size(array);//sizeof(array)/sizeof(array[0]);
  int size = array.size();
  cout<<"  arrayyy     "<<size<<endl;
   findMaximum(array, size, &index,&max);
   array.clear();
  //  delete [] array;
   cout<<" Largest "<<max << "   "<< index<< endl;  
}

// double array[7];
// Entry loop
// sector loop
// J loop
// array[j]=adc[i][j];
// }

// j largest index, + max 
// 				    h[i][index]->max
// 				    del array[]
// 				    } endl i loop
				    
























//void findMaximum(double *num, int size, int *index, double *max)
 void findMaximum( vector<double>num, int size, int *index, double *max)
{
  // int num[5]=
  for (int i=0;i<size;i++)
    {
      if(num[0]<num[i])
	{
	num[0]=num[i];
	*index=i;
	*max=num[0];
	}
      else
	{
	  *max=num[0];
          *index=0;
	}
    }
 
  // return num[0];
}
