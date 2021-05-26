#include <bits/stdc++.h>

using namespace std;
const long double EPS = 1e-7;

class Point {
public:
    long double x, y;
    int num;
    Point() {};

    Point(long double x_, long double y_) {
        x = x_;
        y = y_;
    }

    friend istream &operator>>(istream  &input, Point &p) {
         input >> p.x >> p.y;
         return input;
    }

    friend ostream &operator<<(ostream  &output, Point &p) {
         output  << '(' << p.x  << ',' << p.y << ')';
         return output;
    }

    friend bool operator<(const Point& l, const Point& r)
    {
        if (l.x == r.x) {
            return l.y < r.y;
        } else {
            return l.x < r.x;
        }
    }

    friend bool operator==(const Point& l, const Point& r)
    {
        return l.x == r.x && l.y == l.y;
    }
};

class Face;
class Edge;

class Face {
public:
    Point p;

    Face(Point p_) {
        p = p_;
    }
    vector<Edge*> edges;
};

class Edge {
public:
    Face *f1;
    Face *f2;
    long double a, b, c;
    long double dy, dx;

    Edge(Face *f1_, Face *f2_) {
        f1 = f1_;
        f2 = f2_;
        Point p1 = Point((f1->p.x + f2->p.x)/2., (f1->p.y + f2->p.y)/2.);
        dx = f2->p.x - p1.x;
        dy = f2->p.y - p1.y;

        swap(dx, dy);
        dy = -dy;
        Point p2(p1.x + dx, p1.y + dy);
        a = p1.y - p2.y;
        b = p2.x - p1.x;
        c = -a*p1.x - b * p1.y;
    }

    friend bool isUpper(Edge* e, Point p1, Point p2) {
        long double dx = p2.x - p1.x;
        long double dy = p2.y - p1.y;

        // rotate
        swap(dx, dy);
        dy = -dy;

        p1.x += dx;
        p1.y += dy;

        return e->a*p1.x + e->b*p1.y + e->c < 0;
    }

    Edge *revEdge;
};

bool cmp(const Edge* e1, const Edge* e2)
{
    return atan2(e1->dy, e1->dx) > atan2(e2->dy, e2->dx);
}


Point* intersect(Edge* e1, Edge* e2) {
    if (abs(e1->a*e2->b - e2->a*e1->b) < EPS) {
        return nullptr;
    }

    Point *res = new Point();
    res->x = -(e1->c*e2->b - e2->c*e1->b)/(e1->a*e2->b - e2->a*e1->b);
    res->y = -(e1->a*e2->c - e2->a*e1->c)/(e1->a*e2->b - e2->a*e1->b);

    return res;
}

long double area(Point p1, Point p2, Point p3) {
    return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

long double area(Face* f1, Face* f2, Face* f3) {
    return area(f1->p, f2->p, f3->p);
}

pair<Edge*, Edge*> createEdge(Face *f1, Face *f2)
{
    Edge *e1 = new Edge(f1, f2);
    Edge *e2 = new Edge(f2, f1);
    e1->revEdge = e2;
    e2->revEdge = e1;
    return make_pair(e1, e2);
}

long double dist(Face * f1, Face * f2) {
    return (f1->p.x - f2->p.x) * (f1->p.x - f2->p.x) + (f1->p.y - f2->p.y) * (f1->p.y - f2->p.y);
}

long double dist(Point *p1, Point *p2) {
    return (p1->x - p2->x) * (p1->x - p2->x) + (p1->y - p2->y) * (p1->y - p2->y);
}

pair<int, int> getDownLine(vector<Face*> leftFaces, vector<Face*> rightFaces) {
    int left = 0;
    int right = 0;

    while (area(rightFaces[right], leftFaces[left], leftFaces[(left+1)%leftFaces.size()]) > 0 ||
           area(rightFaces[right], leftFaces[left], leftFaces[(left-1+leftFaces.size())%leftFaces.size()]) > 0 ||
           area(leftFaces[left], rightFaces[right], rightFaces[(right-1+rightFaces.size())%rightFaces.size()]) < 0 ||
           area(leftFaces[left], rightFaces[right], rightFaces[(right+1)%rightFaces.size()]) < 0 ||
           (area(rightFaces[right], leftFaces[left], leftFaces[(left+1)%leftFaces.size()]) == 0 && dist(rightFaces[right], leftFaces[left]) > dist(rightFaces[right], leftFaces[(left+1)%leftFaces.size()])) ||
           (area(leftFaces[left], rightFaces[right], rightFaces[(right-1+rightFaces.size())%rightFaces.size()]) == 0 && dist(leftFaces[left], rightFaces[right]) > dist(leftFaces[left], rightFaces[(right-1+rightFaces.size())%rightFaces.size()])) ||
           (area(rightFaces[right], leftFaces[left], leftFaces[(left-1+leftFaces.size())%leftFaces.size()]) == 0 && dist(rightFaces[right], leftFaces[left]) > dist(rightFaces[right], leftFaces[(left-1+leftFaces.size())%leftFaces.size()])) ||
           (area(leftFaces[left], rightFaces[right], rightFaces[(right+1)%rightFaces.size()]) == 0 && dist(leftFaces[left], rightFaces[right]) > dist(leftFaces[left], rightFaces[(right+1)%rightFaces.size()]))
           ) {
        if (area(rightFaces[right], leftFaces[left], leftFaces[(left-1+leftFaces.size())%leftFaces.size()]) == 0 && dist(rightFaces[right], leftFaces[left]) > dist(rightFaces[right], leftFaces[(left-1+leftFaces.size())%leftFaces.size()])) {
            left--;
            left += leftFaces.size();
            left %= leftFaces.size();
        } else if (area(leftFaces[left], rightFaces[right], rightFaces[(right+1)%rightFaces.size()]) == 0 && dist(leftFaces[left], rightFaces[right]) > dist(leftFaces[left], rightFaces[(right+1)%rightFaces.size()])) {
            right++;
            right %= rightFaces.size();
        } else if (area(rightFaces[right], leftFaces[left], leftFaces[(left+1)%leftFaces.size()]) == 0 && dist(rightFaces[right], leftFaces[left]) > dist(rightFaces[right], leftFaces[(left+1)%leftFaces.size()])) {
            left++;
            left %= leftFaces.size();
        } else if (area(leftFaces[left], rightFaces[right], rightFaces[(right-1+rightFaces.size())%rightFaces.size()]) == 0 && dist(leftFaces[left], rightFaces[right]) > dist(leftFaces[left], rightFaces[(right-1+rightFaces.size())%rightFaces.size()])) {
            right--;
            right += rightFaces.size();
            right %= rightFaces.size();
        } else if (area(rightFaces[right], leftFaces[left], leftFaces[(left+1)%leftFaces.size()]) > 0) {
            left++;
            left %= leftFaces.size();
        } else if (area(rightFaces[right], leftFaces[left], leftFaces[(left-1+leftFaces.size())%leftFaces.size()]) > 0) {
            left--;
            left += leftFaces.size();
            left %= leftFaces.size();
        } else if (area(leftFaces[left], rightFaces[right], rightFaces[(right-1+rightFaces.size())%rightFaces.size()]) < 0) {
            right--;
            right += rightFaces.size();
            right %= rightFaces.size();
        } else if (area(leftFaces[left], rightFaces[right], rightFaces[(right+1)%rightFaces.size()]) < 0) {
            right++;
            right %= rightFaces.size();
        }
    }

    return make_pair(left, right);
}

pair<int, int> getUpLine(vector<Face*> leftFaces, vector<Face*> rightFaces) {
    /*cout << "LEFT\n";
    for (auto i:leftFaces) {
        cout << i->p << ' ';
    }
    cout << '\n';
    cout << "RIGHT\n";
    for (auto i:rightFaces) {
        cout << i->p << ' ';
    }
    cout << '\n';*/

    int left = 0;
    int right = 0;

    while (area(rightFaces[right], leftFaces[left], leftFaces[(left+1)%leftFaces.size()]) < 0 ||
           area(rightFaces[right], leftFaces[left], leftFaces[(left-1+leftFaces.size())%leftFaces.size()]) < 0 ||
           area(leftFaces[left], rightFaces[right], rightFaces[(right-1+rightFaces.size())%rightFaces.size()]) > 0 ||
           area(leftFaces[left], rightFaces[right], rightFaces[(right+1)%rightFaces.size()]) > 0 ||
           (area(rightFaces[right], leftFaces[left], leftFaces[(left-1+leftFaces.size())%leftFaces.size()]) == 0 && dist(rightFaces[right], leftFaces[left]) > dist(rightFaces[right], leftFaces[(left-1+leftFaces.size())%leftFaces.size()])) ||
           (area(leftFaces[left], rightFaces[right], rightFaces[(right+1)%rightFaces.size()]) == 0 && dist(leftFaces[left], rightFaces[right]) > dist(leftFaces[left], rightFaces[(right+1)%rightFaces.size()])) ||
           (area(rightFaces[right], leftFaces[left], leftFaces[(left+1)%leftFaces.size()]) == 0 && dist(rightFaces[right], leftFaces[left]) > dist(rightFaces[right], leftFaces[(left+1)%leftFaces.size()])) ||
           (area(leftFaces[left], rightFaces[right], rightFaces[(right-1+rightFaces.size())%rightFaces.size()]) == 0 && dist(leftFaces[left], rightFaces[right]) > dist(leftFaces[left], rightFaces[(right-1+rightFaces.size())%rightFaces.size()]))
           ) {
        //cout << leftFaces[left]->p << ' ' << rightFaces[right]->p << ' ' << left << ' ' << right << '\n';
        if (area(rightFaces[right], leftFaces[left], leftFaces[(left+1)%leftFaces.size()]) == 0 && dist(rightFaces[right], leftFaces[left]) > dist(rightFaces[right], leftFaces[(left+1)%leftFaces.size()])) {
            left++;
            left %= leftFaces.size();
        } else if (area(leftFaces[left], rightFaces[right], rightFaces[(right-1+rightFaces.size())%rightFaces.size()]) == 0 && dist(leftFaces[left], rightFaces[right]) > dist(leftFaces[left], rightFaces[(right-1+rightFaces.size())%rightFaces.size()])) {
            right--;
            right += rightFaces.size();
            right %= rightFaces.size();
        } else if (area(rightFaces[right], leftFaces[left], leftFaces[(left+1)%leftFaces.size()]) < 0) {
            left++;
            left %= leftFaces.size();
        } else if (area(rightFaces[right], leftFaces[left], leftFaces[(left-1+leftFaces.size())%leftFaces.size()]) == 0 && dist(rightFaces[right], leftFaces[left]) > dist(rightFaces[right], leftFaces[(left-1+leftFaces.size())%leftFaces.size()])) {
            left--;
            left += leftFaces.size();
            left %= leftFaces.size();
        } else if (area(leftFaces[left], rightFaces[right], rightFaces[(right+1)%rightFaces.size()]) == 0 && dist(leftFaces[left], rightFaces[right]) > dist(leftFaces[left], rightFaces[(right+1)%rightFaces.size()])) {
            right++;
            right %= rightFaces.size();
        } else if (area(rightFaces[right], leftFaces[left], leftFaces[(left-1+leftFaces.size())%leftFaces.size()]) < 0) {
            left--;
            left += leftFaces.size();
            left %= leftFaces.size();
        } else if (area(leftFaces[left], rightFaces[right], rightFaces[(right-1+rightFaces.size())%rightFaces.size()]) > 0) {
            right--;
            right += rightFaces.size();
            right %= rightFaces.size();
        } else if (area(leftFaces[left], rightFaces[right], rightFaces[(right+1)%rightFaces.size()]) > 0) {
            right++;
            right %= rightFaces.size();
        }
    }
    //cout << leftFaces[left]->p << ' ' << rightFaces[right]->p << ' ' << left << ' ' << right << '\n';
    //cout << area(leftFaces[left], rightFaces[right], rightFaces[(right+1+rightFaces.size())%rightFaces.size()]) << ' ' << rightFaces[(right+1+rightFaces.size())%rightFaces.size()]->p << '\n';
    //cout << "FOUND UP\n";
    char c;
    //cin >> c;
    return make_pair(left, right);
}

bool onSameSide(Edge* e, Point p1, Point p2) {
    return (e->a*p1.x + e->b*p1.y + e->c >= 0) == (e->a*p2.x + e->b*p2.y + e->c >= 0) ||
    abs(e->a*p2.x + e->b*p2.y + e->c) < EPS ||
    abs(e->a*p1.x + e->b*p1.y + e->c) < EPS;
}

Point* intersect(Face* f1, Face* f2, int edgePos) {
    pair<Edge*, Edge*> edges = createEdge(f1, f2);
    Edge* e = edges.first;
    Point *p = intersect(f1->edges[edgePos], e);
    if (p == nullptr) {
        return nullptr;
    }

    for (auto e:f1->edges) {
        if (!onSameSide(e, f1->p, *p)) {
            return nullptr;
        }
    }

    return p;
}

void printEdges(Face *f) {
    return;
    cout << "------------------------------\n";
    cout << "EDGES FOR " << f->p << '\n';
    for (auto e:f->edges) {
        cout << e->f2->p << '\n';
    }
    cout << "------------------------------\n";
}

void fixEdges(Face* f) {
    sort(f->edges.begin(), f->edges.end(), cmp);

    for (int i = 0; i + 1 < f->edges.size(); i++) {
        if (area(f->edges[i]->f2, f, f->edges[i+1]->f2) <= 0) {
            vector<Edge*> edges;
            for (int j = i + 1; j < f->edges.size(); j++) {
                edges.push_back(f->edges[j]);
            }
            for (int j = 0; j <= i; j++) {
                edges.push_back(f->edges[j]);
            }
            f->edges = edges;
            return;
        }
    }
}

void splitChain(vector<Face*> leftFaces, vector<Face*> rightFaces, pair<int, int> up, pair<int, int> down) {
    Face *left = leftFaces[up.first];
    Face *right = rightFaces[up.second];

    vector<Edge*> addLeft, addRight;
    int leftEdge = left->edges.size() - 1;
    int startLeftEdge = leftEdge;
    int rightEdge = 0;
    int startRightEdge = rightEdge;
    bool isStartLeft = true;
    bool isStartRight = true;
    Point *prevLeft = nullptr;
    Point *prevRight = nullptr;

    while (!(left == leftFaces[down.first] && right == rightFaces[down.second])) {
        //cout << left->p << ' ' << right->p << '\n';
        Point *leftIntersect = nullptr;
        int cnt = 0;
        while (leftIntersect == nullptr && cnt < left->edges.size()) {
            cnt++;
            if (!isStartLeft && leftEdge == startLeftEdge) {
                break;
            }
            leftIntersect = intersect(left, right, leftEdge);
            if (area(right, left, left->edges[leftEdge]->f2 ) < 0) {
                leftIntersect = nullptr;
            }

            /*
            if (leftIntersect != nullptr && area(left->p, *leftIntersect, left->edges[leftEdge]->f2->p)  0) {
                leftIntersect = nullptr;
            }*/
            if (leftEdge == startLeftEdge && prevLeft != nullptr && leftIntersect != nullptr && dist(prevLeft, leftIntersect) < EPS) {
                leftIntersect = nullptr;
            }
            if (addLeft.size() != 0 && leftIntersect != nullptr) {
                for (auto edge: addLeft) {
                    if (((edge->a * left->p.x + edge->b * left->p.y + edge->c) >= 0) != ((edge->a * leftIntersect->x + edge->b * leftIntersect->y + edge->c) >= 0) && abs(edge->a * leftIntersect->x + edge->b * leftIntersect->y + edge->c) > EPS) {
                        //cout << "BBBBBBBBBBBBB\n";
                        leftIntersect = nullptr;
                        break;
                    }
                }
            }
            if (leftIntersect != nullptr) {
                break;
            }
            leftEdge--;
            leftEdge += left->edges.size();
            leftEdge %= left->edges.size();
        }
        Point *rightIntersect = nullptr;
        cnt = 0;
        while (rightIntersect == nullptr  && cnt < right->edges.size()) {
            //cout << "checking " << right->edges[rightEdge]->f2->p << ' ' << rightEdge << '\n';
            cnt++;
            if (!isStartRight && rightEdge == startRightEdge) {
                break;
            }
            rightIntersect = intersect(right, left, rightEdge);
            if (area(left, right, right->edges[rightEdge]->f2 ) > 0) {
                rightIntersect = nullptr;
            }
            //cout << rightIntersect << '\n';
            /*
            if (rightIntersect != nullptr && area(right->p, *rightIntersect, right->edges[rightEdge]->f2->p) < 0) {
                rightIntersect = nullptr;
            }*/

            if (rightEdge == startRightEdge && prevRight != nullptr && rightIntersect != nullptr && dist(prevRight, rightIntersect) < EPS) {
                // cout << *prevRight << ' ' << *rightIntersect << ' ' << dist(prevRight, rightIntersect) << '\n';
                //cout << "WWWWWW\n";
                rightIntersect = nullptr;
            }
            if (addRight.size() != 0 && rightIntersect != nullptr) {
                for (auto edge: addRight) {
                    if (((edge->a * right->p.x + edge->b * right->p.y + edge->c) >= 0) != ((edge->a * rightIntersect->x + edge->b * rightIntersect->y + edge->c) >= 0) && abs(edge->a * rightIntersect->x + edge->b * rightIntersect->y + edge->c) > EPS) {
                        // cout << (edge->a * rightIntersect->x + edge->b * rightIntersect->y + edge->c) << '\n';
                        //cout << "AAAAAAAAAAAAAA\n";
                        rightIntersect = nullptr;
                        break;
                    }
                }
            }
            if (rightIntersect != nullptr) {
                break;
            }
            rightEdge++;
            rightEdge %= right->edges.size();
        }

        //cout << "LEFT EDGE " <<  left->edges[leftEdge]->f2->p << '\n';
        //cout << "RIGHT EDGE " <<  right->edges[rightEdge]->f2->p << '\n';
        //cout << "INTERSECTIONS " << leftIntersect << ' ' << rightIntersect << '\n';
        //cout << *leftIntersect << ' ' << *rightIntersect << '\n';

        pair<Edge*, Edge*> edges = createEdge(left, right);
        addLeft.push_back(edges.first);
        addRight.push_back(edges.second);
        //cout << "AAAA\n";
        //cout << (rightIntersect == nullptr || (leftIntersect != nullptr && isUpper(edges.first, *leftIntersect, *rightIntersect))) << '\n';
        if (rightIntersect == nullptr || (leftIntersect != nullptr && isUpper(edges.first, *leftIntersect, *rightIntersect))) {
            prevLeft = leftIntersect;
            //cout << "BBBB\n";
            vector<Edge*> newEdges;
            for (int j = addLeft.size()-1; j >= 0; j--) {
                newEdges.push_back(addLeft[j]);
            }
            for (int j = leftEdge; true; j=(j+1)%left->edges.size()) {
                newEdges.push_back(left->edges[j]);
                if (j == startLeftEdge) {
                    break;
                }
            }
            Edge *revEdge = left->edges[leftEdge]->revEdge;
            left->edges = newEdges;
            fixEdges(left);
            printEdges(left);
            //cout << "EDGES FOR " << left->p << '\n';
            for (auto j:left->edges) {
                //cout << j->f2->p << '\n';
            }

            addLeft.clear();

            left = revEdge->f1;
            leftEdge = 0;
            while (left->edges[leftEdge] != revEdge) {
                leftEdge++;
            }
            startLeftEdge = leftEdge;
            leftEdge--;
            leftEdge += left->edges.size();
            leftEdge %= left->edges.size();
            rightEdge = startRightEdge;
            rightEdge++;
            rightEdge %= right->edges.size();

            isStartRight = true;
            isStartLeft = false;
        } else {
            //cout << "BBBB\n";
            prevRight = rightIntersect;
            vector<Edge*> newEdges;
            for (int j = startRightEdge; true; j=(j+1+right->edges.size())%right->edges.size()) {
                newEdges.push_back(right->edges[j]);
                if (j == rightEdge) {
                    break;
                }
            }
            //cout << "BBBB\n";
            for (int j = 0; j < addRight.size(); j++) {
                newEdges.push_back(addRight[j]);
            }
            Edge *revEdge = right->edges[rightEdge]->revEdge;
            //cout << "EDGES FOR " << right->p << '\n';
            for (auto j:right->edges) {
                //cout << j->f2->p << '\n';
            }
            //cout << "@@@\n";
            printEdges(right);
            right->edges = newEdges;
            fixEdges(right);
            printEdges(right);

            //cout << "EDGES FOR " << right->p << '\n';
            for (auto j:right->edges) {
                //cout << j->f2->p << '\n';
            }
            addRight.clear();
            //cout << "BBBB\n";

            right = revEdge->f1;
            rightEdge = 0;
            while (right->edges[rightEdge] != revEdge) {
                //cout << right->edges[rightEdge]->f2->p << '\n';
                rightEdge++;
                if (rightEdge == right->edges.size()) {
                    cout << "STOP";
                    exit(1);
                }
            }
            //cout << "BBBB\n";
            startRightEdge = rightEdge;
            rightEdge++;
            rightEdge %= right->edges.size();
            leftEdge = startLeftEdge;
            leftEdge--;
            leftEdge += left->edges.size();
            leftEdge %= left->edges.size();
            isStartLeft = true;
            isStartRight = false;
        }
        char c;
        //cin >> c;
    }
    //cout << left->p << ' ' << right->p << '\n';
    //cout << "END\n";
    fflush(stdout);
    pair<Edge*, Edge*> edges = createEdge(left, right);
    addLeft.push_back(edges.first);
    addRight.push_back(edges.second);

    vector<Edge*> newEdges;
    for (int j = addLeft.size()-1; j >= 0; j--) {
        newEdges.push_back(addLeft[j]);
    }
    for (int j = 0; j<=startLeftEdge; j++) {
        newEdges.push_back(left->edges[j]);
    }
    left->edges = newEdges;
    fixEdges(left);
    printEdges(left);

    //cout << "END2\n";
    fflush(stdout);

    newEdges.clear();
    for (int j = startRightEdge; j<right->edges.size(); j++) {
        newEdges.push_back(right->edges[j]);
    }
    for (int j = 0; j < addRight.size(); j++) {
        newEdges.push_back(addRight[j]);
    }
    printEdges(right);
    right->edges = newEdges;
    fixEdges(right);
    printEdges(right);
    //cout << "EDGES FOR " << right->p << '\n';
    for (auto j:right->edges) {
        //cout << j->f2->p << '\n';
    }

    fflush(stdout);
}

vector<Face*> mergeFaces(vector<Face*> leftFaces, vector<Face*> rightFaces) {
    pair<int, int> up = getUpLine(leftFaces, rightFaces);
    pair<int, int> down = getDownLine(leftFaces, rightFaces);

    splitChain(leftFaces, rightFaces, up, down);

    if (up.first == down.first && up.second == down.second) {
        vector<Face*> res;

        int cnt = leftFaces.size();
        int cur = down.first;
        while (cnt--) {
            res.push_back(leftFaces[cur]);
            cur++;
            cur %= leftFaces.size();
        }
        reverse(res.begin(), res.end());

        cnt = rightFaces.size();
        cur = up.second;
        while (cnt--) {
            res.push_back(rightFaces[cur]);
            cur++;
            cur %= rightFaces.size();
        }

        return res;
    }

    vector<Face*> res;
    int cur = down.first;
    while (cur != up.first) {
        res.push_back(leftFaces[cur]);
        cur++;
        cur %= leftFaces.size();
    }
    res.push_back(leftFaces[cur]);

    cur = up.second;
    while (cur != down.second) {
        res.push_back(rightFaces[cur]);
        cur++;
        cur %= rightFaces.size();
    }
    res.push_back(rightFaces[cur]);

    return res;
}
Face* faces[100000];

vector<Face*> solve(vector<Point> vec) {
    if (vec.size() == 1) {
        vector<Face*> res;
        Face *f = new Face(vec[0]);
        faces[vec[0].num] = f;
        res.push_back(f);
        return res;
    }

    if (vec.size() == 2) {
        vector<Face*> res;
        Face *f1 = new Face(vec[0]);
        Face *f2 = new Face(vec[1]);
        faces[vec[0].num] = f1;
        faces[vec[1].num] = f2;
        pair<Edge*, Edge*> edges = createEdge(f1, f2);
        f1->edges.push_back(edges.first);
        f2->edges.push_back(edges.second);
        res.push_back(f1);
        res.push_back(f2);
        return res;
    }

    if (vec.size() == 3) {
        vector<Face*> res;
        Face *f1 = new Face(vec[0]);
        Face *f2 = new Face(vec[1]);
        Face *f3 = new Face(vec[2]);

        faces[vec[0].num] = f1;
        faces[vec[1].num] = f2;
        faces[vec[2].num] = f3;
        if (abs(area(vec[0], vec[1], vec[2])) < EPS) {
            pair<Edge*, Edge*> edges = createEdge(f1, f2);
            f1->edges.push_back(edges.first);
            f2->edges.push_back(edges.second);

            edges = createEdge(f3, f2);
            f3->edges.push_back(edges.first);
            f2->edges.push_back(edges.second);
        } else {
            if (area(f1, f2, f3) < 0) {
                swap(f2, f3);
            }

            pair<Edge*, Edge*> edges = createEdge(f1, f2);
            f1->edges.push_back(edges.first);
            f2->edges.push_back(edges.second);

            edges = createEdge(f3, f2);
            f3->edges.push_back(edges.first);
            f2->edges.push_back(edges.second);

            edges = createEdge(f3, f1);
            f3->edges.push_back(edges.first);
            f1->edges.push_back(edges.second);
            reverse(f1->edges.begin(), f1->edges.end());
        }

        res.push_back(f1);
        res.push_back(f3);
        res.push_back(f2);
        return res;
    }

    vector<Point> left, right;

    for (int i = 0; i < vec.size(); i++) {
        if (i < vec.size()/2) {
            left.push_back(vec[i]);
        } else {
            right.push_back(vec[i]);
        }
    }

    vector<Face*> leftFaces = solve(left);
    vector<Face*> rightFaces = solve(right);
    return mergeFaces(leftFaces, rightFaces);
}

void debug() {
/*
(9075,1377)
(19430,11480)
(16430,10392)
(31051,21968)
*/
    vector<Point> vec = {Point(9075, 1377), Point(19430, 11480), Point(15, 15), Point(10, 12), Point(9, 11), Point(3, 3), Point(6, 2), Point(16, 10), Point(11, 5), Point(13, 19)};
    sort(vec.begin(), vec.end());
    solve(vec);
}

void test() {
    srand(0);

    int n;
    cin >> n;
    int cnt = 0;
    int md = 100000;

    for (int i = 1; i <= 1000; i++) {
        cnt++;
        cout << cnt << '\n';
        vector<Point> vec;
        set<int> setx, sety;
        for (int i = 0; i < n; i++) {
            Point p;
            p.x = rand() % md;
            p.y = rand() % md;
            while (!(setx.find(p.x) == setx.end() && sety.find(p.y) == sety.end())) {
                p.x = rand() % md;
                p.y = rand() % md;
            }
            setx.insert(p.x);
            sety.insert(p.y);
            //cout << p << '\n';
            vec.push_back(p);
        }
        sort(vec.begin(), vec.end());
        vector<Face*> res = solve(vec);
    }
}

void rot(Point *p) {
    long double x = p->x*cos(1) - p->y*sin(1);
    long double y = p->x*sin(1) + p->y*cos(1);

    p->x = x;
    p->y = y;
}

int main()
{
    //test();
    //debug();
    freopen("C:\\Users\\38067\\Desktop\\ogkg\\bin\\Debug\\input.txt", "r", stdin);
    freopen("C:\\Users\\38067\\Desktop\\ogkg\\bin\\Debug\\output.txt", "w", stdout);
    //srand(time(0));
    vector<Point> vec, vecc;
    int n;
    cin >> n;
    map<pair<int, int>, int> st;
    for (int i = 1; i <= n; i++) {
        Point p;
        cin >> p.x;
        cin >> p.y;
        if (st[make_pair(p.x, p.y)]) {
            continue;
        }
        st[make_pair(p.x, p.y)]++;
        p.num = i;
        vecc.push_back(p);
        rot(&p);
        vec.push_back(p);
    }

    sort(vec.begin(), vec.end());
    solve(vec);

    for (int i = 1; i <= n; i++) {
        long double best = 0;
        for (auto e:faces[i]->edges) {
            if (best == 0 || dist(&vecc[i-1], &vecc[e->f2->p.num-1]) < dist(&vecc[i-1], &vecc[best-1])) {
                best = e->f2->p.num;
            }
        }
        if (st[make_pair(vecc[i-1].x, vecc[i-1].y)] != 1) {
            cout << i-1 << '\n';
        } else {
            cout << best-1 << '\n';
        }
    }
}

/*
(9,8)
(7,6)
(0,2)
(6,1)
(4,3)
(2,5)
(3,7)
(1,4)
(8,9)
*/
