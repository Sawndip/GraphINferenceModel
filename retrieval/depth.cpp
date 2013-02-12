#include <iostream>
#include <vector>

using namespace std;

void pv(vector<int> v)
{
    for(vector<int>::iterator it = v.begin(); it != v.end(); it++)
    {
	cout << *it << " ";
    }
    cout << endl;
}

int main() {
    int x = 10; 
    int y = 10;
    
    int vx = --x;
    int vy = y--;
    
    vector<int> v;
    
    for(int i=0; i < 10; i++) 
	v.push_back(i);
    
    pv(v);
    v.erase(v.begin());
    pv(v);
   
    
}