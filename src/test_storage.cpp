#include <bits/stdc++.h>

using namespace std;

unsigned int *arr;


void add_to_arr(int x, int j, int beta){

    unsigned temp = j*beta;
    unsigned idx = (temp) >> 5;
    unsigned rem=temp&(0x0000001F); 

    unsigned add = rem+beta;

    if(add > 32){
        arr[idx] = arr[idx] & ( ((1 << rem) - 1) << (32-rem) );
        arr[idx] = arr[idx] | (x >> (add-32));
        arr[idx+1] = (x << (64-add)) | ( arr[idx+1] & ( (1<<64-add) - 1)   )  ;
    }
    else{
        unsigned mask_left,mask_right;
        mask_left = ( (1<<rem) - 1) << (32-rem);
        mask_right = (1<<(32-add)) -1 ;
        arr[idx] = (arr[idx] & mask_left) | (arr[idx] & mask_right);
        arr[idx] = arr[idx] | (x<<(32-add)) ;
    } 
    
}

int extract_from_arr(int j,int beta){

    unsigned temp = j*beta;
    unsigned idx = (temp) >> 5;
    unsigned rem=temp&(0x0000001F);
    
    int res=0;
    
    unsigned add = rem + beta;

    if(add > 32){
        unsigned  mask=0;
        mask = (1 << (32-rem)) - 1;
        int lastXbits = arr[idx] & mask;
        res = (lastXbits << (add-32)) | (arr[idx+1] >> (64-add));       
    }
    else{
        res = (arr[idx]<<rem)>>(32-beta);
    }
    return res;
}


int *p;
int main(int argc, char * argv[]){

    int max_bits = atoi(argv[1]);
    int size_p =188380 * 50;

    int max_size = (int)pow(2,max_bits);


    srand (time(NULL));

    p = new int[size_p];

    int NInts = int(max_bits*size_p/32) + 1;
    arr = new unsigned int[NInts];

   
    for(int i=0;i<size_p;i++){
        p[i] = rand() % (max_size);
        add_to_arr(p[i],i+1,max_bits);
    }

    // cout<<"P array"<<endl;
    // for(int i=0;i<size_p;i++)
    //     cout<<p[i] << "  ";


    // cout<<"\nP array in binary"<<endl;
    // for(int i=0;i<size_p;i++)
    //     cout<<bitset<max_bits>(p[i]);

    // cout<<"\nP array in binary space efficient"<<endl;
    // for(int i=0;i<NInts;i++)
    //     cout<<bitset<32>(arr[i]);

    // cout <<"\nOriginal p[i] and stored p[i]"<<endl;

    for(int i=0;i<size_p;i+=17){
        p[i] = rand() % (max_size);
        add_to_arr(p[i],i+1,max_bits);
    }

    int max_diff=0;
    for(int i=0;i<size_p;i++){
        int temp = extract_from_arr(i+1,max_bits);
        max_diff = max(max_diff,abs(p[i]-temp));
    }
    cout<<"Error: " << max_diff<<endl;
}