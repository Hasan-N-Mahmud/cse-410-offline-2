
#include<windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <bits/stdc++.h>
#include <stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include <sstream>

using namespace std;
const double PI = 3.14159;
struct point
{
	double x,y,z,w;
	point(){
        x=y=z=0;
        w=1;
        }
	point(double a,double b,double c){
	    x=a;
	    y=b;
	    z=c;
	    w=1;
	}
	point(double a,double b,double c,double d){
	    x=a;
	    y=b;
	    z=c;
	    w=d;
	}

	point operator+ (point p) {
	    point res;
	    res.x = x + p.x;
	    res.y = y + p.y;
	    res.z = z + p.z;
	    return res; }
    void print(){
        cout<<x<<"\t"<<y<<"\t"<<z<<"\t"<<w<<endl;
    }
    void normalize(){
        if(w!=1){
            x /= w;
            y /= w;
            z /= w;
            w /= w;
        }
    }
    void precision(){
        x =roundf(x * 10000) / 10000.0 ;
        y =roundf(y * 10000) / 10000.0;
       z = roundf(z * 10000) / 10000.0;
    }

};

struct triangle
{
	point a,b,c;
	triange(){

        }
	triangle(point x,point y,point z){
	    a=x;
	    b=y;
	    c=z;
	}
friend ostream &operator<<( ostream &output, const triangle &t ) {
         output <<  t.a.x << " " << t.a.y << " "<<t.a.z<<"\n"<< t.b.x << " " << t.b.y << " "<<t.b.z<<"\n"<<  t.c.x << " " << t.c.y << " "<<t.c.z<<"\n";
         return output;
      }
    void print(){
        cout<<"triangle printing"<<endl;
        a.print();
        b.print();
        c.print();
        cout<<endl;
    }
    void precision(){
        a.precision();
        b.precision();
        c.precision();
        }
};

struct matrix{
    int row,col;
    double * mat;
    matrix(int a,int b){
    row = a;
    col = b;
    mat = new double[row*col];
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
                {
                    mat[i*col+j] =0;
                }
    }
    }

    matrix(int a,int b,string str){
    row = a;
    col = b;
    mat = new double[row*col];
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
                {
                    if(i == j)
                        mat[i*col+j] =1;
                    else
                    mat[i*col+j] =0;
                }
    }
    }
    void print(){
    for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
                {
                    cout << mat[i*col+j] << " ";
                }
            cout << endl;
    }
}

matrix operator*(matrix m2){
    matrix res(row,col);
    int i, j, k;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            for (k = 0; k < col; k++)
                res.mat[i*col+j] += mat[i*col+k] * m2.mat[k*col+j];
        }
    }
    return res;
}

point operator*(point p){
    matrix res(4,1);
    double arr[4];
    arr[0]=p.x;
    arr[1]=p.y;
    arr[2]=p.z;
    arr[3]=p.w;
    int i, j;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
                res.mat[i] += mat[i*col+j] * arr[j];
        }
    }
    point r(res.mat[0],res.mat[1],res.mat[2],res.mat[3]);
    r.normalize();
    return r;
}
bool operator==(matrix m2){
    int i,j;
        for (i = 0; i < row; i++) {
            for (j = 0; j < col; j++) {
                if(mat[i*col+j] != m2.mat[i*col+j])
                    return false;
        }
    }
    return true;
}
bool is_identity_matrix(){
    int i,j;
        for (i = 0; i < row; i++) {
            for (j = 0; j < col; j++) {
                if((i==j && mat[i*col+j] !=1) || (i != j && mat[i*col+j] !=0))
                    return false;
        }
    }
    return true;
}
};

point point_init(string s)
{
    stringstream ss(s);
    double word;
    vector<double> values;
    while (ss >> word) {
        values.push_back(word);
    }
   point temp(values[0],values[1],values[2]);
   return temp;
}

struct Vector{
    double x,y,z;
    Vector(){};
    Vector(double a,double b, double c){
        x=a;
        y=b;
        z=c;
    }
    Vector operator+ (Vector p) {
	    Vector res;
	    res.x = x + p.x;
	    res.y = y + p.y;
	    res.z = z + p.z;
	    return res; }
Vector operator- (Vector p) {
	    Vector res;
	    res.x = x - p.x;
	    res.y = y - p.y;
	    res.z = z - p.z;
	    return res; }
    Vector operator*(double k) {
	    Vector res;
	    res.x = x * k;
	    res.y = y * k;
	    res.z = z * k;
	    return res; }
	    void normalize(){
            double value = sqrt(x*x + y*y + z*z);
            x /=  value;
            y /=  value;
            z /=  value;
}
	    void print(){
	    cout<<"Vector printing..."<<endl;
	    cout<<"x: "<<x<<" y: "<<y<<" z: "<<z<<endl;
	    }
	    void scrutinize(){
        if(x == -0){
            x+=1;
            x-=1;
	    }
	    if(y == -0){
            y+=1;
            y-=1;
	    }
	    if(z == -0){
            z+=1;
            z-=1;
	    }
	    }
};
double dotProduct(Vector vect_A, Vector vect_B)
{
    double product = 0;
    product += vect_A.x * vect_B.x;
    product += vect_A.y * vect_B.y;
    product += vect_A.z * vect_B.z;
    product+=1;
    product -=1;
    return product;
}

Vector crossProduct(Vector vect_A, Vector vect_B)
  {
     Vector cross_P;
     cross_P.x = vect_A.y * vect_B.z - vect_A.z * vect_B.y;
     cross_P.y = vect_A.z * vect_B.x - vect_A.x * vect_B.z;
     cross_P.z = vect_A.x * vect_B.y - vect_A.y * vect_B.x;
     cross_P.scrutinize();
     return cross_P;
}

Vector vector_init(string s)
{
    stringstream ss(s);
    double word;
    vector<double> values;
    while (ss >> word) {
        values.push_back(word);
    }
   Vector temp(values[0],values[1],values[2]);
   return temp;
}

Vector Rod(Vector x,Vector a,double angle){
    Vector res;
    double cosValue,sinValue;
    cosValue = cos(angle * PI /180.0);
    sinValue = sin(angle * PI /180.0);
    cout << fixed << setprecision(4) << cosValue<<" "<<sinValue <<endl;
    res = x * cosValue + a * ((1 - cosValue)*dotProduct(a,x))  ;
    res = res +  crossProduct(a,x) * sinValue;
    return res;
}
Vector eye,look,up;
double fovY,aspectRatio,near_,far_;
int main()
{
matrix current_trans(4,4,"identity");
//current_trans.print();
stack<matrix> trans_stack;
stack<matrix> ins_history;
vector<triangle> triangles;
bool push_flag = false;
string myText;
ifstream MyReadFile("scene.txt");
int lineCount=1,i,j,k=1;
ins_history.push(current_trans);
while (getline (MyReadFile, myText)) {

    switch(lineCount){
    case 1:
        eye = vector_init(myText);
        //eye.print();
        break;
    case 2:
        look = vector_init(myText);
       // look.print();
        break;
    case 3:
        up = vector_init(myText);
        //up.print();
        break;
    case 4:
        stringstream ss(myText);
        double word;
        vector<double> values;
        while (ss >> word) {
            values.push_back(word);
        }
        fovY = values[0];
        aspectRatio = values[1];
        near_ = values[2];
        far_ = values[3];
        break;
    }
    if(lineCount > 4){
        if(myText == "triangle"){
            vector<point> temp;
            for(int i=0;i<3;i++){
            getline (MyReadFile, myText);
                temp.push_back(point_init(myText));
            }
            triangle a(temp[0],temp[1],temp[2]);
            //cout<<"Current Trans Matrix"<<endl;
            //current_trans.print();
            a.a = current_trans * a.a;
            a.b = current_trans * a.b;
            a.c =  current_trans * a.c;
            a.precision();
            triangles.push_back(a);
            //cout<<"Triangle: "<<k<<endl;
            k++;
            //a.print();
        }
    else if(myText == "scale"){
        cout<<"scale\n";
        getline (MyReadFile, myText);
        matrix t(4,4);
        stringstream ss(myText);
        double sValue;
        i=0,j=0;
        while(ss >> sValue){
            t.mat[i*4+j ] = sValue;
            i++;
            j++;
        }
        t.mat[i*4+j]=1;
        //t.print();
        if(current_trans.is_identity_matrix()){
            current_trans = t;
            //cout<<"init Current trans\n";
            //current_trans.print();
        }
        else{
            current_trans = current_trans * t;
             //cout<<"init Current trans\n";
            //current_trans.print();
        }
        ins_history.push(current_trans);
        if(push_flag){
            trans_stack.push(current_trans);
            //cout<<"top\n";
            //trans_stack.top().print();
            push_flag = false;
        }
    }
        else if(myText == "translate"){
        cout<<"translate\n";
        getline (MyReadFile, myText);
        matrix t(4,4,"identity");
        stringstream ss(myText);
        double tValue;
        i=0;
        while(ss >> tValue){
            t.mat[i*4+3]=tValue;
            i++;
        }
        //t.print();
        if(current_trans.is_identity_matrix()){
            current_trans = t;
             //cout<<"init Current trans\n";
             //current_trans.print();
        }
        else{
            current_trans = current_trans * t;
             //cout<<" Current trans\n";
             //current_trans.print();
        }
        ins_history.push(current_trans);
        //cout<<"top history\n";
        //ins_history.top().print();
        if(push_flag){
            trans_stack.push(current_trans);
            //cout<<"top\n";
            //trans_stack.top().print();
            push_flag = false;
        }
    }
    else if(myText == "rotate"){
        cout<<"Rotate"<<endl;
        getline (MyReadFile, myText);
        matrix t(4,4,"identity");
        stringstream ss(myText);
        double angle;
        Vector a,c1,c2,c3;
        ss >> angle;
        ss >> a.x;
        ss >> a.y;
        ss >> a.z;
        a.normalize();
        c1 = Rod(Vector(1,0,0),a,angle);
        c2 = Rod(Vector(0,1,0),a,angle);
        c3 = Rod(Vector(0,0,1),a,angle);
    t.mat[0] = c1.x;
    t.mat[1] = c2.x;
    t.mat[2] = c3.x;
    t.mat[3] = 0;
    t.mat[4] = c1.y;
    t.mat[5] = c2.y;
    t.mat[6] = c3.y;
    t.mat[7] = 0;
    t.mat[8] = c1.z;
    t.mat[9] = c2.z;
    t.mat[10] = c3.z;
    t.mat[11] = 0;
    //cout<<"Rotation Matrix\n";
    //t.print();
    //cout<<endl;
     if(current_trans.is_identity_matrix()){
            current_trans = t;
             //cout<<"init Current trans\n";
             //current_trans.print();
        }
        else{
            current_trans = current_trans * t;
             //cout<<" Current trans\n";
             //current_trans.print();
        }
        ins_history.push(current_trans);
        //cout<<"top history\n";
        //ins_history.top().print();
        if(push_flag){
            trans_stack.push(current_trans);
            //cout<<"top\n";
            //trans_stack.top().print();
            push_flag = false;
        }
    }
    else if(myText == "push"){
        cout<<"push\n";
        push_flag = true;
    }
    else if(myText == "pop"){
        cout<<"pop\n";
        while(!(ins_history.top() == trans_stack.top())){
            ins_history.pop();
        }
        //cout<<"check\n";
       // trans_stack.top().print();
        ins_history.pop();
        trans_stack.pop();
        current_trans = ins_history.top();
        //cout<<"After pop\n";
        //current_trans.print();
    }
    else if(myText == "end"){
        break;
    }

    }
    lineCount++;
}

MyReadFile.close();
string fName = "stage1.txt";
ofstream MyFile(fName);

for(i=0;i<triangles.size();i++){
    triangle t = triangles.at(i);
    //triangles.pop_back();
    MyFile<<t<<endl;
}
MyFile.close();
//View Transformation

Vector l,r,u;
l = look - eye;
l.normalize();
r = crossProduct(l,up);
r.normalize();
u = crossProduct(r,l);

matrix T_matrix(4,4,"identity");
T_matrix.mat[3]= -eye.x;
T_matrix.mat[7]= -eye.y;
T_matrix.mat[11]= -eye.z;
matrix R_matrix(4,4,"identity");
R_matrix.mat[0] = r.x;
R_matrix.mat[1] = r.y;
R_matrix.mat[2] = r.z;
R_matrix.mat[4] = u.x;
R_matrix.mat[5] = u.y;
R_matrix.mat[6] = u.z;
R_matrix.mat[8] = -l.x;
R_matrix.mat[9] = -l.y;
R_matrix.mat[10] = -l.z;
matrix V_matrix = R_matrix * T_matrix;
MyFile.open("stage2.txt");
for(i=0;i<triangles.size();i++){
    triangle t = triangles.at(i);
    t.print();
    t.a = V_matrix * t.a;
    t.b = V_matrix * t.b;
    t.c = V_matrix * t.c;
    t.print();
    cout<<"Vcetor\n";
    triangles.insert(triangles.begin()+i,t);
    triangles.erase(triangles.begin()+i+1);
    triangles.at(i).print();
    MyFile<<t<<endl;
}
MyFile.close();
return 0;
}
