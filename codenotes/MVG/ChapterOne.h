//ChapterOne.h
#ifndef CODENOTES_MVG_CHAPTERONE_H
#define CODENOTES_MVG_CHAPTERONE_H

namespace codenotes {

namespace mvg {

/**********************
 * 第一章 2D射影几何和变换
 * 
 * 第一章主要介绍平面的射影变换几何，这些变换模拟的是透视摄像机对平面摄像时所产生  \
 * 的几何失真；
 * 
 * 本章的目的是从透视图像中恢复仿射性质（平行线等）和度量性质（线的夹角等）；
 **********************/


/**************************
 * 1-1
 * 点、直线、二次曲线的表示形式；
 * 射影变换下，几何实体的映射；
 * ************************
 * IR^2：平面的矢量空间；
 * IP^2：IR^3 - (0, 0, 0)， 即射影空间，通过射影点所有直线的集合，  \
 * 减去(0, 0, 0)是因为0矢量不与任何直线相对应;
 **************************/
 


/****************************
 * 直线表示
 * 直线方程：ax+by+c=0
 * 所以直线可以使用矢量(a,b,c)来表示；
 * 直线乘上一个任意非0的缩放因子，并不发生改变，所以平面的直线的自由度是2；
 * 
 * 点表示
 * 由直线方程推出：(x,y,1)(a,b,c)^T = (x,y,1)I = 0
 * 可以看出2D点和线在表达形式上非常相似，我们同样可以用(x1,x2,x3)来表示一个点，与直线类似；
 * 
 * 性质
 * 1. 点是否在直线上的判断：x^T * I = 0， 即两个向量的点积；
 * 2. 两直线的交点：l和l'的叉乘，叉乘表示与两个直线相垂直的直线，也就是与两个直线的点积为0,可以将
 *    此直线看作一个点；
 * 3. 过两点的直线：两点的叉积；
 *****************************/

//其实直线和点完全可以用一个类来进行表达，但是为了符合常识，这里还是分别进行实现；
//直线实现
//我认为直线的参数应该限制为浮点型，所以这里只实现浮点型参数的直线；
template<typename T> class Point2;

class Line2f {
public:
  Line2f():representation({0, 0}) { }

  Line2f(float v1, float v2, float v3)
    :representation({v1/v3, v2/v3}) { }

public:
  template<typename T>
  inline bool throughPoint(Point2<T> p) const {
      return (representation[0] * p.x() + representation[1] * p.y() + 1) == 0;
  }

  inline Point2f intersectionWithLine(Line2f l) const {
      float k = representation[0] * l.intercept() - representation[2] * l.slope();
      return Point2f((representation[1] - l.intercept()) / k,  \
                     (l.slope() - representation[0]) / k);  //返回值优化
  }

  inline float slope() const {
      return representation[0];
  }

  inline float intercept() const {
      return representation[1];
  }

private:
  float representation[2];  //矢量表示(representation[0], representation[1], 1)
};

//点实现
template<typename T>
class Point2 {
public:
  Point2():representation({0, 0}) { }

  Point2(T v1, T v2)
    :representation({v1, v2}) { }

public:
  inline bool inLine(Line2f l) const {
      return (representation[0] * l.slope() + representation[1] * l.intercept() + 1) == 0;
  }

  inline Line2f lineWithPoint(Point2<T> p) const {
      return Line2f(representation[1] - p.y(), p.x() - representation[0],  \
                    representation[0] * p.y() - representation[2] * p.x());  //返回值优化
  }

  inline T x() const {
      return representation[0];
  }

  inline T y() const {
      return representation[1];
  }

private:
  T representation[2];  //矢量表示(representation[0], representation[1], 1)
};

typedef Point2<float>    Point2f;
typedef Point2<int>    Point2i;
typedef Point2<unsigned int>    Point2ui;

/*****************************
 * 理想点和无穷远线
 * 平面上两个不平行的线必相交于一点；但是平行的线并不相交，或者说交于无穷远，可以用(x1,x2,0)表示;
 * 这些无穷点称为理想点，我们发现这些理想点都在一条线上，(0,0,1);
 * 任意两条平行的线都在无穷远相交同一点；
 * 对偶定理：2D射影几何的任何定理都有一个对应的对偶定理，可以通过互换点和线的作用导出；
 * 
 * 从理想点和无穷远线角度考虑，上面Point2和Line2f的实现只考虑了欧式有限的点，而没有无穷远的表达；
 ******************************/

/*******************************
 * 二次曲线
 * 非齐次方程：ax^2+bxy+cy^2+dx+ey+f=0
 * x->x1/x3  y->x2/x3进行齐次化:ax1^2+bx1x2+cx2^2+dx1x3+ex2x3+fx3^2=0
 * ====>>矩阵形式： x^T * C * x = 0
 *        C = [a    b/2    d/2
 *             b/2    c    e/2
 *             d/2    e/2    f]
 * 二次曲线具有5个自由度
 * 
 * 二次曲线C上x的切线I由 I = Cx确定
 * 证明：I^T * x = x^T * C * x = 0  ->   I=Cx  (C是对称矩阵)
 *      如果I和C只交于一点，则可证明I是C的曲线；
 *      假设I与C交于另一点y，则y^T * C * y =0； x^T * C * y = I^T * y =0;
 *      则(x+ay)^T * C * (x+ay) = 0
 *      这说明x,y连线上都在C，如果C没有退化成直线，则说明I只与C交于一点；
 * 
 * 对偶二次曲线
 * 5个点可以确定二次曲线，那么5个切线确定的二次曲线为对偶形式，I^T*C‘*I=0
 * C’为C的伴随矩阵，非奇异矩阵C'=C.inv
 * *******************************/

/*********************************
 * 射影变换
 * 射影变换也就是我们常说的单应性变换；
 * 透视变换，是一种特殊的射影变换，常见于相机成像过程中；
 * 透视变换和射影变换主要不同是：透视变换中两个平面对应的点的连线都经过同一个点，而射影变换则不是；
 * 
 * 点的映射：x'=Hx
 * 线的映射：I'=inv(H)^T*I
 * 
 * 性质：射影变换保证了变换前同一条直线上的点变换后仍然在同一条线上；
 * 证明：直线I上的点x，有I^T * x = 0,
 *      则I^T * inv(H) * H * x = 0;
 *      所以变换后点都落在inv(H)^T * I直线上；
 * 
 * 
 * 
 * 平行线射影变换后不再保持平行
 * I1 = （a, b, c1)  I2= (a, b, c2)
 * x1 = cross_product(I1, I2) = (b*c2-b*c1, c1*a - a*c2, 0)  //I1和I2的交点；
 * 
 * I1' = H * I1 = (h11*a+h12*b+h13*c1, h21*a+h22*b+h23*c1, h31*a+h32*b+h33*c1), 
 * I2' = H * I2 = (h11*a+h12*b+h13*c2, h21*a+h22*b+h23*c2, h31*a+h32*b+h33*c2), 
 * 
 * x2 = cross_product(I1', I2')
 * x2(2) = (h11*a+h12*b+h13*c1)*(h21*a+h22*b+h23*c2)-
 *         (h21*a+h22*b+h23*c1)*(h11*a+h12*b+h13*c2)
 *       =  (X+h13*c1)*(Y+h23*c2)-(Y+h23*c1)*(X+h13*c2)
 *       = XY+X*h23*c2+Y*h13*c1+h13*h23*c1*c2 -(XY +Y*h13*c2+X*h23*c1+h13*h23*c1*c2)
 *       = X*h23*c2+Y*h13*c1-Y*h13*c2-X*h23*c1
 * 
 * x2(2)!=0,所以映射后不再保持平行
 * 
 *
 *  
 * 应用：消除透视失真
 * 如果可以计算出来H矩阵，再利用inv(H)则可以实现失真(平行保持)的校正；
 * 理论上通过4对两个平面对应的坐标，实现H的计算,但是透视映射只有6个自由度，所以需要3个对应点即可；
 * 但是需要人工选择，是否有什么方法让它自动校正呢？
 *************************************/



/*************************************
 * 几种常见的平面变换
 * 欧式变换：只进行旋转和平移，长度和面积都保持不变；
 * 相似变换：在欧式变换的基础上，增加了缩放系数，长度比，夹角和虚圆点都保持不变；
 * 仿射变换：角度会被改变，平行，面积比和平行线段的长度比不变；（
 * 射影变换：只有直线的性质被保留下来；
 * 
 * 仿射变换通过SVD进行分解，可以分解为两个旋转以及非均匀的拉伸变换；
 * 
 * 
 * 射影变换与仿射变换相比，变换矩阵最后一行多了(v1,v2,1)，这也让射影变换可以将无穷远点映射到有限点的能力；
 * 
 * 根据自己的需求，可以选择矫正到相应的变换即可。
 ***************************************/



/****************************************
 * 1-7 从图像恢复仿射和度量性质
 * 
 * 如果我们要确定长度比，只需要矫正到相似变换即可，只需要四个自由度，可以由无穷远线和线上的两个虚圆点确定；
 * x' = Hp*x
 * x''=Hp'*x'
 * x'' = Ha*x
 * Ha=Hp*Hp'
 * (0,0,1)=Hp'*(l1,l2,l3)^T
 * 
 * 
 * 
 * 在射影变换H下，无穷远线为不动直线的充要条件是H是仿射变换；
 * 
 * 由图像恢复仿射性质
 * 
 * 经过射影变换后两条平行的直线相交于一点，如果有两个这样的点，就可以确定一条由无穷远直线映射的有限直线；
 * 为了得到仿射性质，只需要将这条有限直线映射到无穷远线就可以了；
 * H=Ha*[1 0 0
 *       0 1 0
 *       l1 l2 l3]
 * Ha是射影变换，H是仿射变换
 * 
 * 如果知道长度比的两个线段，则可以由直线上的距离比确定无穷远点；
 * 
 * 
 * 在射影变换H下，虚圆点I和J为不动点的充要条件是H是相似变换；
 * 这里需要先将射影变换矫正到仿射变换，再将仿射变换矫正到相似变换；
 * 
 * 由虚圆点可以确定与虚圆点对应的对偶二次曲线，C* = I*J^T+J*I^T;
 * 对偶二次曲线在射影变换H下不变的充要条件是H是相似变换；
 * 二次曲线确认后，同样可以计算出欧式交，长度比；
 * 
 * 如何确定虚原点？
 * 如果确定了椭圆，通过椭圆与无穷远线的交点就是虚圆点；
 * 利用两组正交线进行矫正；
 ********************************************/


}  //mvg

}  //codenotes

#endif
