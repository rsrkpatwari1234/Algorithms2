/*
NAME : RADHIKA PATWARI
ROLL NO. : 18CS10062
Lab Assignment no. 2 of Algorithms 2 for Semester 5

RUN COMMAND : gcc asgn2_18CS10062.c -o a -lm
			  ./a < input.txt
*/

//including the C header files 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <math.h>

#define PI 3.14159265358979323846			//predefined value for PI

//structure for Edge
typedef struct POINT{
	double x;
	double y;
}POINT;

//structure for Line Segment
typedef struct LINESEG{
	double startX, startY;
	double endX, endY;
}LINESEG;

//structure for Arc
typedef struct ARC{
	double centerX, centerY;
	double fromAngle, toAngle;
}ARC;

// Merge sort implementation
// Merges two subarrays of S[].  
// First subarray is arr[l..m]  
// Second subarray is arr[m+1..r]  
void merge(POINT *S, int l, int m, int r)  
{  
    int i, j, k;  
    int n1 = m - l + 1;  
    int n2 = r - m;  
    POINT L[n1], R[n2];  			// create temp arrays
 
    for (i = 0; i < n1; i++)  		// Copy data to temp array L[]
        L[i] = S[l + i];  
    for (j = 0; j < n2; j++)  		// Copy data to temp array R[]
        R[j] = S[m + 1 + j];  
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray  
    j = 0; // Initial index of second subarray  
    k = l; // Initial index of merged subarray  
    while (i < n1 && j < n2) {  
        if (L[i].x < R[j].x) {  		//sorted in ascending order of x coordinates
            S[k] = L[i];  
            i++;  
        }  
        else if (L[i].x > R[j].x){  
            S[k] = R[j];  
            j++;  
        }
        else{							//for the case of equal x coordinates
        	if (L[i].y > R[j].y){       // use descending order of y coordinates
	            S[k] = L[i];  
	            i++;  
	        } 
	        else{  
	            S[k] = R[j];  
	            j++;  
	        }
        }  
        k++;  
    }  
    while (i < n1) {  	//Copy the remaining elements of L[], if there are any 
        S[k] = L[i];  
        i++;  
        k++;  
    }  
    while (j < n2) {    //Copy the remaining elements of R[], if there are any
        S[k] = R[j];  
        j++;  
        k++;  
    }  
}  
  
//l is for left index and r is right index of the sub-array of S to be sorted
void mergeSort(POINT *S, int l, int r)  
{  
    if (l < r) {  
        int m = l + (r - l) / 2;  // Same as (l+r)/2, but avoids overflow for large l and h
        mergeSort(S, l, m);  	  // Sort first half of array of points
        mergeSort(S, m + 1, r);   // Sort second half of array of points
        merge(S, l, m, r);  	  // Merging the two halves of array of points	
    }  
}  

// Printing the points present in pts array of points
void PrintPoints(int n, POINT *pts){
	for(int i = 0; i < n; i++){
		printf("   %.15lf	%.15lf\n", pts[i].x, pts[i].y);
	}
}

// Computing the slope for 2 points p1 and p2 as slope=(p2.y-p1.y)/(p2.x-p1.x)
double slopeComputation(POINT first, POINT second){
	double slope = (second.y - first.y) / (1.0 * (second.x - first.x));
	return slope;
}

// Computing orientation of p3 wrt line segment p1p2
// returning true if popping needs to be done (p3 is collinear or lies on left of segment p1p2)
// returning false if no popping is required (p3 lies in right of the line segment p1p2)
bool orientation(POINT p1, POINT p2, POINT p3){

	double value1 = (1.0) * (p2.y - p1.y)*(p3.x - p2.x);
	double value2 = (1.0) * (p2.x - p1.x)*(p3.y - p2.y); 

	return (value1 > value2) ? false : true;
}

// Graham Scan algorithm for computing Convex Hull in O(nlogn) time and O(n) space
// S : sorted order of points in ascending x coordinate and descending y coordinate if required
// flag : 0 for finding UPPER HULL
// flag : 1 for finding LOWER HULL
// returning # of points on the respective HULL
int CH(POINT *S, int n, int flag, POINT *H){

	if(n <= 2)					//sanity check for possiblity of convex hull
		return 0;				//Minimum possible Convex hull is of size 3

	int start_point, next_point, i, size;
	double slope, final_slope;

	if(flag){					//for LOWER HULL case,reverse S 
		int j = n-1;
		POINT temp;
		for(i = 0; i < j; i++, j--){
			temp = S[i];
			S[i] = S[j];
			S[j] = temp;
		}
	}

	start_point = 0;			//take the leftmost and rightmost point for UPPER and LOWER HULL respectively
	next_point = 1;				//initialise 2nd point
	
	// points with same x coordinate do not actually exist
	while(S[next_point].x == S[start_point].x)   
		next_point++;

	//finding next point with least angle from the vertical ray of start_point in clockwise direction
	final_slope = slopeComputation(S[start_point], S[next_point]);

	for(i = next_point+1; i < n; i++){
		slope = slopeComputation(S[start_point], S[i]);		//slope calculation
		if(slope >= final_slope){
			final_slope = slope;
			next_point = i;  
		}
	}

	H[0] = S[start_point];			// setting first point in the HULL array
	H[1] = S[next_point];  			// setting second point in the HULL array

	size = 1;

	// if the point lies on right of line segment,insert it in the HULL array
	// if the point is collinear with the line segment,remove topmost ele from HULL array
	// if the point lies on left of line segment,remove topmost ele from HULL array
	for(i = next_point+1; i<n; i++){
		while(orientation(H[size-1], H[size], S[i])){  //collinear or left side case
			size--;				
		}
		size++;
		H[size] = S[i];
	}
	if(flag){						//reversing back S for LOWER HULL case
		int j = n-1;
		POINT temp;
		for(i = 0; i < j; i++, j--){
			temp = S[i];
			S[i] = S[j];
			S[j] = temp;
		}
	}
	return size+1;
}

// computing the  segment and arcs of the boundary using the obtained Convex Hulls
// applying geometry for computing angles of arcs and start and end coordinates of tangents
void contzone(POINT *UH, int u, POINT *LH, int l, double radius, LINESEG *T, ARC *A){

	int i;									
	double slopePrev, slopeNext;	// slopePrev : slope between previous and current points 
									// slopeNext : slope between current and next points
	//for upper hull
	for(i = 0; i < u; i++){
		if(i == 0){					//leftmost point of the UPPER HULL
			slopeNext = slopeComputation(UH[0], UH[1]);
			A[0].centerX = UH[0].x;
			A[0].centerY = UH[0].y;
			A[0].fromAngle = PI;
			A[0].toAngle = PI/2.0 + atan(slopeNext);
			T[0].startX = UH[0].x + (radius * cos(A[0].toAngle));
			T[0].startY = UH[0].y + (radius * sin(A[0].toAngle));
			slopePrev = slopeNext;
		}	
		else if(i == u-1){			//rightmost point of the UPPER HULL
			A[u-1].centerX = UH[u-1].x;
			A[u-1].centerY = UH[u-1].y;
			A[u-1].fromAngle = PI/2.0 + atan(slopePrev);
			A[u-1].toAngle = 0.0;
			T[i-1].endX = UH[u-1].x + (radius * cos(A[u-1].fromAngle));
			T[i-1].endY = UH[u-1].y + (radius * sin(A[u-1].fromAngle));
		}
		else{
			slopeNext = slopeComputation(UH[i], UH[i+1]);
			A[i].centerX = UH[i].x;
			A[i].centerY = UH[i].y;
			A[i].fromAngle = PI/2.0 + atan(slopePrev);
			A[i].toAngle = PI/2.0 + atan(slopeNext);
			T[i-1].endX = UH[i].x + (radius * cos(A[i].fromAngle));
			T[i-1].endY = UH[i].y + (radius * sin(A[i].fromAngle));
			T[i].startX = UH[i].x + (radius * cos(A[i].toAngle));
			T[i].startY = UH[i].y + (radius * sin(A[i].toAngle));
			slopePrev = slopeNext;
		}
	}

	//for lower hull
	for(i = 0; i < l; i++){
		if(i == 0){					//rightmost point of the CONVEX HULL
			slopeNext = slopeComputation(LH[0], LH[1]);
			A[u].centerX = LH[0].x;
			A[u].centerY = LH[0].y;
			A[u].fromAngle = 0.0;
			A[u].toAngle = atan(slopeNext)-(PI/2.0);
			T[u-1].startX = LH[0].x - (radius * cos(PI+A[u].toAngle));
			T[u-1].startY = LH[0].y - (radius * sin(PI+A[u].toAngle));
			slopePrev = slopeNext;
		}	
		else if(i == l-1){			//leftmost point of the CONVEX HULL
			A[u+l-1].centerX = LH[l-1].x;
			A[u+l-1].centerY = LH[l-1].y;
			A[u+l-1].fromAngle = atan(slopePrev)-(PI/2.0);
			A[u+l-1].toAngle = -1 * PI;
			T[u-2+l-1].endX = LH[l-1].x - (radius * cos(PI+A[u+l-1].fromAngle));
			T[u-2+l-1].endY = LH[l-1].y - (radius * sin(PI+A[u+l-1].fromAngle));
		}
		else{
			slopeNext = slopeComputation(LH[i], LH[i+1]);
			A[u+i].centerX = LH[i].x;
			A[u+i].centerY = LH[i].y;
			A[u+i].fromAngle = atan(slopePrev)-(PI/2.0);
			A[u+i].toAngle = atan(slopeNext)-(PI/2.0);
			T[u-2+i].endX = LH[i].x - (radius * cos(PI+A[u+i].fromAngle));
			T[u-2+i].endY = LH[i].y - (radius * sin(PI+A[u+i].fromAngle));
			T[u-1+i].startX = LH[i].x - (radius * cos(PI+A[u+i].toAngle));
			T[u-1+i].startY = LH[i].y - (radius * sin(PI+A[u+i].toAngle));
			slopePrev = slopeNext;
		}
	}
}

// Printing the segments and arcs for each of the upper and lower hulls
void printcontzone(int u, int l, LINESEG *T, ARC *A){
	int i;
	printf("\n+++ The containment zone\n");
	printf("--- Upper section\n");
	for(i=0;i<u-1;i++){
		printf("    Arc      : (%.15lf,%.15lf) From %.15lf to %.15lf \n",A[i].centerX, A[i].centerY, A[i].fromAngle, A[i].toAngle);
		printf("    Tangent  : From (%.15lf,%.15lf) to (%.15lf,%.15lf) \n",T[i].startX, T[i].startY, T[i].endX, T[i].endY);
	}
	printf("    Arc      : (%.15lf,%.15lf) From %.15lf to %.15lf \n",A[u-1].centerX, A[u-1].centerY, A[u-1].fromAngle, A[u-1].toAngle);
	printf("--- Lower section\n");
	for(i=0;i<l-1;i++){
		printf("    Arc      : (%.15lf,%.15lf) From %.15lf to %.15lf \n",A[u+i].centerX, A[u+i].centerY, A[u+i].fromAngle, A[u+i].toAngle);
		printf("    Tangent  : From (%.15lf,%.15lf) to (%.15lf,%.15lf) \n",T[u-1+i].startX, T[u-1+i].startY, T[u-1+i].endX, T[u-1+i].endY);
	}
	printf("    Arc      : (%.15lf,%.15lf) From %.15lf to %.15lf \n",A[u+l-1].centerX, A[u+l-1].centerY, A[u+l-1].fromAngle, A[u+l-1].toAngle);

}

// Main() function
int main(){

	//input from terminal
	int n;
	double radius;
	scanf("%d", &n);
	scanf("%lf", &radius);

	POINT *S;
	S = (POINT *)malloc(n * sizeof(POINT));	//dynamic memory allocation for n points

	double x,y;
	int i;

	for(i = 0; i < n; i++){
		scanf("%lf", &x);
		scanf("%lf", &y);
		S[i].x = x;
		S[i].y = y;
	}

	printf("%d\n", n);					// printing # of points
	printf("%.15lf\n", radius);			// printing radius of each point
	for(i = 0; i < n; i++){							
		printf("%.15lf	%.15lf\n", S[i].x, S[i].y);	// printing each point
	}

	mergeSort(S, 0, n-1);				// sorting points in ascending x coordinate 
										// and descending y coordinate for same x case
	printf("\n+++ Circles after sorting\n");
	PrintPoints(n, S);

	// computing UPPER HULL
	POINT *upperH;
	upperH = (POINT *)malloc(n * sizeof(POINT));
	int upperN;
	upperN = CH(S, n, 0, upperH);
	printf("\n+++ Upper hull\n");
	PrintPoints(upperN, upperH);

	// computing LOWER HULL
	POINT *lowerH;
	lowerH = (POINT *)malloc(n * sizeof(POINT));
	int lowerN;
	lowerN = CH(S, n, 1, lowerH);
	printf("\n+++ Lower hull\n");
	PrintPoints(lowerN, lowerH);

	// computing LINE SEGMENTS and ARCS of the boundary of contaminated region
	LINESEG* Tangents;
	ARC *Arcs;
	Tangents = (LINESEG *)malloc(n * sizeof(LINESEG));
	Arcs = (ARC *)malloc((upperN+lowerN) * sizeof(ARC));
	contzone(upperH, upperN, lowerH, lowerN, radius, Tangents, Arcs);
	printcontzone(upperN, lowerN, Tangents, Arcs);

	// releasing the dynamically allocated space
	free(S);					
	free(upperH);
	free(lowerH);
	free(Tangents);
	free(Arcs);

	return 0;	
}