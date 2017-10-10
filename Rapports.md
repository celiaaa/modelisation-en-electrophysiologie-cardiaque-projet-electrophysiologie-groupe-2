cd# Rapports hebdomadaires -- à faire en fin de séance
## Mardi 3/10/2017

Première partie : recherche de resources, documentation  

On remarque que:
(\partial_t u - \partial_x(a(x) \partial u) = f(t,x)) est une équation du type Sturm-Liouville (qu'on notera SL):   

Les équation de type SL:(\partial_x(a(x) \partial u) + q(x)*u = \lamda*r(x)*u  

http://www.bibnum.education.fr/mathematiques/analyse/le-probleme-de-sturmliouville  
https://maths.ucd.ie/~ekashdan/page2/SLP.pdf  
https://www.math.univ-toulouse.fr/~rondep/CoursTD/polyic2.pdf  
http://www.math.ubc.ca/~israel/m316/nhslp.pdf  


Méthodes des différences finis pour la résolution des équations SL:  

http://www.sciencedirect.com/science/article/pii/089812219390323N  
http://www.ams.org/journals/mcom/1969-23-108/S0025-5718-1969-0258291-7/S0025-5718-1969-0258291-7.pdf  



1.si a(x) est égal à une constante k on aura une equation de la chaleur non homogène et on a comme condition au bord:
u(t,0)=u⁰
u(t,1)=u¹
pour tout t appartenant à (0,l1)
u(0,x)=u⁰(x) pour tout x appartenant à (0,l2)

Conditions aux limites problème de Sturm-Trouville du type : 
c1y(a) + c2y'(a) = 0 ; c1y(b) + c2y'(b) = 0
on prendra c1=d1=0

si a(x) dépend de x on aura comme conditions limites
pour x appartenant à l'interval [a,b] avec a<b on a ;

partial_xu =h*u=0 si x=a;
partial_xu + H*u=0 si x=b;   et on a h,H appartienne à ]0,+infini[
u(0,x)=g(x)
il faudra aussi que a(x)>0;

## Mardi 17/10/2017
## Mardi 14/11/2017
## Mardi 21/11/2017
## Mardi 28/12/2017
