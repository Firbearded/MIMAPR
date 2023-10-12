#include <vector>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <fstream>

using namespace std;


class Point
{
public:
	double x, y;
	Point(double x_, double y_) : x(x_), y(y_) {}
	Point() : x(0), y(0) {}
};


class EdgeCondition  // Класс граничного условия
{
	int type;
	double T;
public:
	int get_type() { return type; }
	int get_T() { return T; }
	EdgeCondition(int type_, double T_ = 0) : type(type_), T(T_) {}
};


class Form  // Класс геометрической формы
{
private:
	EdgeCondition cond;
public:
	virtual bool is_inside(double x, double y, bool is_add) = 0;
	virtual double get_dist_x(double x, double y) = 0;
	virtual double get_dist_y(double x, double y) = 0;
	virtual char get_form() = 0;

	virtual Point get_p1() = 0;
	virtual Point get_p2() = 0;

	Form(EdgeCondition c) : cond(c) {}
	EdgeCondition get_cond() { return cond; }
};


class Rectangle : public Form  // Класс квадрата
{
	Point p1, p2;
public:
	Rectangle(double x1, double y1, double x2, double y2, EdgeCondition c) : Form(c)
	{
		p1.x = min(x1, x2);
		p1.y = min(y1, y2);
		p2.x = max(x1, x2);
		p2.y = max(y1, y2);
	}
	bool is_inside(double x, double y, bool is_add)
	{
		if (is_add)
			return x > p1.x && x < p2.x&& y > p1.y && y < p2.y;
		else
			return x >= p1.x && x <= p2.x&& y >= p1.y && y <= p2.y;
	}
	double get_dist_x(double x, double y)  // Надо перевернуть обратно, возможно у тоже и просмотреть их использование
	{
		if (abs(x - p1.x) < abs(x - p2.x))
			return p1.x - x;
		else
			return p2.x - x;
	}
	double get_dist_y(double x, double y)
	{
		if (abs(y - p1.y) < abs(y - p2.y))
		{
			//cout << p1.y - y << endl;
			return p1.y - y;
		}
		else
		{
			//cout << p2.y - y << endl;
			return p2.y - y;
		}
	}
	char get_form() { return 'r'; }
	Point get_p1() { return p1; }
	Point get_p2() { return p2; }
};


class Circle : public Form  // Класс круга
{
	double r;
	Point o;
public:
	Circle(double x, double y, double r_, EdgeCondition c) : Form(c), o(Point(x, y)), r(r_) {}
	bool is_inside(double x, double y, bool is_add)
	{
		if (is_add)
			return (o.x - x) * (o.x - x) + (o.y - y) * (o.y - y) < r * r;
		else
			return (o.x - x) * (o.x - x) + (o.y - y) * (o.y - y) <= r * r;
	}
	double get_dist_x(double x, double y)
	{
		if (abs(-x + o.x - sqrt(abs(r * r - (y - o.y) * (y - o.y)))) < abs(-x + o.x + sqrt(abs(r * r - (y - o.y) * (y - o.y)))))
			return -x + o.x - sqrt(abs(r * r - (y - o.y) * (y - o.y)));
		else
			return -x + o.x + sqrt(abs(r * r - (y - o.y) * (y - o.y)));
	}
	double get_dist_y(double x, double y)
	{
		if (abs(o.y - sqrt(abs(r * r - (x - o.x) * (x - o.x))) - y) < abs(o.y + sqrt(abs(r * r - (x - o.x) * (x - o.x))) - y))
		{
			//cout << o.y - sqrt(abs(r * r - (x - o.x) * (x - o.x))) - y << endl;
			return o.y - sqrt(abs(r * r - (x - o.x) * (x - o.x))) - y;
		}
		else
		{
			//cout << o.y + sqrt(abs(r * r - (x - o.x) * (x - o.x))) - y << endl;
			return o.y + sqrt(abs(r * r - (x - o.x) * (x - o.x))) - y;
		}
	}
	char get_form() { return 'c'; }
	Point get_p1() { return Point(o.x - r, o.y - r); }
	Point get_p2() { return Point(o.x + r, o.y + r); }
};


class Node;


class EdgeNode : public Point  // Класс граничного узла
{
	double *T;
	double mu;
	double h, h_t;
	Node* parent;
	EdgeCondition cond;
public:
	EdgeNode(Node *p_, double x, double y, EdgeCondition c, double mu_, int sim_t, double h_, double h_t_, double T_ = 0) :
		parent(p_), cond(c), Point(x, y), mu(mu_), h(h_), h_t(h_t_)
	{
		T = new double[sim_t + 1];
		T[0] = T_;
	}
	void get_edge_dif(double* difw, double* difr)
	{
		switch (cond.get_type())
		{
		case 3:
			*difr = -(cond.get_T() * 2 * h_t / (h * h * mu * (1 + mu)));
			*difw = 0;
			break;
		case 4:
			*difr = cond.get_T() * 2 * h_t / (h * (1 + mu));
			*difw = 2 * h_t / (h * h * mu * (1 + mu));
			break;
		case 5:
			*difr = 0;
			*difw = 2 * h_t / (h * h * mu * (1 + mu) * (1 + h * mu));
			break;
		}
	}
	void calculate(int t);
	double get_mu() { return mu; }
	double get_T(int t) { return T[t]; }
	EdgeCondition get_cond() { return cond; }
};


class Node : public Point  // Класс не граничного узла
{
	bool active;
	Form* parent;
	double* T;
	EdgeNode* edge_x;
	EdgeNode* edge_y;
public:
	
	Node(double x, double y, double T_, Form *p, bool active_, int sim_t) : Point(x, y), parent(p), edge_x(NULL), edge_y(NULL), active(active_)
	{
		T = new double[sim_t + 1];
		T[0] = T_;
	}
	bool is_active() { return active; }
	Form* get_parent() { return parent; }
	EdgeNode* get_edge_x() { return edge_x; };
	EdgeNode* get_edge_y() { return edge_y; };
	double get_T(int t) { return T[t]; }
	void set_active(bool a) { active = a; }
	void set_parent(Form* p) { parent = p; }
	void set_edge_x(EdgeNode* ex) { edge_x = ex; }
	void set_edge_y(EdgeNode* ey) { edge_y = ey; }
	void set_T(int t, double T_) { T[t] = T_; }
};

void EdgeNode::calculate(int t)
{
	switch (cond.get_type())
	{
	case 3:
		T[t] = cond.get_T();
		break;
	case 4:
		T[t] = parent->get_T(t) - mu * h * cond.get_T();
		break;
	case 5:
		T[t] = (mu * parent->get_T(t) / (mu * h + 1));
		break;
	default:
		cout << "\n\nERROR!\n\n";
		exit(0);
		break;
	}
}


class Object  // Класс пластины
{
private:
	vector<Form*> forms;
	int n, m;
	int sim_t;
	double real_n, real_m;
	double h, h_t;
	bool** bin_mask;
	Node*** nodes;
	vector<EdgeNode*> edges;
public:
	Object(double real_n_, double real_m_, double h_, double sim_t_, double h_t_) :
		n(real_n_ / h_ + 1), m(real_m_ / h_ + 1), real_n(real_n_), real_m(real_m_), h(h_), sim_t(sim_t_ / h_t_), h_t(h_t_)
	{
		bin_mask = new bool* [m];
		nodes = new Node * *[m];
		for (int i = 0; i < m; i++)
		{
			bin_mask[i] = new bool[n];
			nodes[i] = new Node * [n];
			for (int j = 0; j < n; j++)
			{
				bin_mask[i][j] = 0;
				nodes[i][j] = NULL;
			}
		}
	}
	~Object()
	{
		for (int i = 0; i < m; i++)
		{
			delete[] bin_mask[i];
			for (int j = 0; j < n; j++)
				delete nodes[i][j];
			delete[] nodes[i];
		}
		delete[] bin_mask;
		delete[] nodes;
		for (int i = 0; i < edges.size(); i++)
			delete edges[i];
	}

	void add_form(Form* f, double NU, bool is_add)  // Добавление не граничных узлов в пластину по геометрической форме
	{
		forms.push_back(f);
		for (int i = (f->get_p1()).y / h; i <= (f->get_p2()).y / h + 1; i += 1)
			for (int j = (f->get_p1()).x / h; j <= (f->get_p2()).x / h + 1; j += 1)
			{
				if (i >= 0 && i * h <= real_m && j >= 0 && j * h <= real_n && f->is_inside(j * h, i * h, is_add))
				{
					if (nodes[i][j])
						delete nodes[i][j];
					nodes[i][j] = new Node(j, i, NU, f, is_add, sim_t);
					bin_mask[i][j] = is_add;
				}
			}
	}

	void compile_edges()  // Добавление граничных узлов
	{
		Form* p;
		double dist, mu;
		int di[4] = { -1, 0, 0, 1 }, dj[4] = { 0, 1, -1, 0 };

		for (int i = 0; i < m; i++)  // Удаление мелких мю, превращение их в дырки с граничным условием как у детали, из которой удалили
			for (int j = 0; j < n; j++)
			{
				if (!nodes[i][j] || !nodes[i][j]->is_active())
				{
					for (int k = 0; k < 4; k++)
					{
						if (i + di[k] < 0 || i + di[k] >= m || j + dj[k] < 0 || j + dj[k] >= n)
							continue;
						if (bin_mask[i + di[k]][j + dj[k]])
						{
							if (nodes[i][j] && !nodes[i][j]->is_active())
								p = nodes[i][j]->get_parent();
							else
								p = nodes[i + di[k]][j + dj[k]]->get_parent();
							dist = min(abs(p->get_dist_x((j + dj[k]) * h, (i + di[k]) * h)), abs(p->get_dist_y((j + dj[k]) * h, (i + di[k]) * h)));
							if (dist <= 1 / h + 0.0000001 && dist >= 0.0000001)
							{
								//elements[i + di[k]][j + dj[k]]->set_parent(p);
								nodes[i + di[k]][j + dj[k]]->set_active(0);
							}
						}
					}
					
				}
			}
		//*/
		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
			{
				if (!nodes[i][j] || !nodes[i][j]->is_active())
				{
					if (i != 0 && bin_mask[i - 1][j] && nodes[i - 1][j]->is_active())
					{
						if (nodes[i][j] && !nodes[i][j]->is_active())
							p = nodes[i][j]->get_parent();
						else
							p = nodes[i - 1][j]->get_parent();
						if (bin_mask[i][j])
							dist = h;
						else
							dist = p->get_dist_y(j * h, (i - 1) * h);
						
						edges.push_back(new EdgeNode(nodes[i - 1][j], j * h, (i - 1) * h + dist, p->get_cond(), abs(dist / h), sim_t, h, h_t));
						nodes[i - 1][j]->set_edge_y(edges[edges.size() - 1]);
					}
					if (i != m - 1 && bin_mask[i + 1][j] && nodes[i + 1][j]->is_active())
					{
						if (nodes[i][j] && !nodes[i][j]->is_active())
							p = nodes[i][j]->get_parent();
						else
							p = nodes[i + 1][j]->get_parent();
						if (bin_mask[i][j])
							dist = h;
						else
							dist = p->get_dist_y(j * h, (i + 1) * h);
						
						if (j != 0 && bin_mask[i][j - 1] && nodes[i][j - 1]->get_edge_x() && (nodes[i][j - 1]->get_edge_x())->x == j * h
							&& (nodes[i][j - 1]->get_edge_x())->y == (i + 1) * h + dist)
							nodes[i + 1][j]->set_edge_y(nodes[i][j - 1]->get_edge_x());
						else if (j != n - 1 && bin_mask[i][j + 1] && nodes[i][j + 1]->get_edge_x() && (nodes[i][j + 1]->get_edge_x())->x == j * h
							&& (nodes[i][j + 1]->get_edge_x())->y == (i + 1) * h + dist)
							nodes[i + 1][j]->set_edge_y(nodes[i][j + 1]->get_edge_x());
						else
						{
							edges.push_back(new EdgeNode(nodes[i + 1][j], j * h, (i + 1) * h + dist, p->get_cond(), abs(dist / h), sim_t, h, h_t));
							nodes[i + 1][j]->set_edge_y(edges[edges.size() - 1]);
						}
					}
					if (j != 0 && bin_mask[i][j - 1] && nodes[i][j - 1]->is_active())
					{
						if (nodes[i][j] && !nodes[i][j]->is_active())
							p = nodes[i][j]->get_parent();
						else
							p = nodes[i][j - 1]->get_parent();
						if (bin_mask[i][j])
							dist = h;
						else
							dist = p->get_dist_x((j - 1) * h, i * h);
						
						if (i != 0 && bin_mask[i - 1][j] && nodes[i - 1][j]->get_edge_y() && (nodes[i - 1][j]->get_edge_y())->y == i * h
							&& (nodes[i - 1][j]->get_edge_y())->x == (j - 1) * h + dist)
							nodes[i][j - 1]->set_edge_x(nodes[i - 1][j]->get_edge_y());
						else
						{
							edges.push_back(new EdgeNode(nodes[i][j - 1], (j - 1) * h + dist, i * h, p->get_cond(), abs(dist / h), sim_t, h, h_t));
							nodes[i][j - 1]->set_edge_x(edges[edges.size() - 1]);
						}
					}
					if (j != n - 1 && bin_mask[i][j + 1] && nodes[i][j + 1]->is_active())
					{
						if (nodes[i][j] && !nodes[i][j]->is_active())
							p = nodes[i][j]->get_parent();
						else
							p = nodes[i][j + 1]->get_parent();
						if (bin_mask[i][j])
							dist = h;
						else
							dist = p->get_dist_x((j + 1) * h, i * h);
						
						if (i != 0 && bin_mask[i - 1][j] && nodes[i - 1][j]->get_edge_y() && (nodes[i - 1][j]->get_edge_y())->y == i * h
							&& (nodes[i - 1][j]->get_edge_y())->x == (j + 1) * h + dist)
							nodes[i][j + 1]->set_edge_x(nodes[i - 1][j]->get_edge_y());
						else
						{
							edges.push_back(new EdgeNode(nodes[i][j + 1], (j + 1) * h + dist, i * h, p->get_cond(), abs(dist / h), sim_t, h, h_t));
							nodes[i][j + 1]->set_edge_x(edges[edges.size() - 1]);
						}
					}
				}
			}
	}
	void calculate_explicit()  // Вычисление значений температур в узлах явным методом
	{
		for (int i = 0; i < edges.size(); i++)
			edges[i]->calculate(0);
		for (int t = 1; t <= sim_t; t++)
		{
			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++)
					if (nodes[i][j] && nodes[i][j]->is_active())
					{
						nodes[i][j]->set_T(t, 0);
						EdgeNode* ex = nodes[i][j]->get_edge_x();
						EdgeNode* ey = nodes[i][j]->get_edge_y();
						double sum = nodes[i][j]->get_T(t - 1);
						if (ex && ex->x < (nodes[i][j]->x) * h)
							sum += h_t * 2 * ((ex->get_mu() * nodes[i][j + 1]->get_T(t - 1) - (ex->get_mu() + 1) * nodes[i][j]->get_T(t - 1) + ex->get_T(t - 1))
								/ (h * h * ex->get_mu() * (ex->get_mu() + 1)));
						else if (ex)
							sum += h_t * 2 * ((ex->get_mu() * nodes[i][j - 1]->get_T(t - 1) - (ex->get_mu() + 1) * nodes[i][j]->get_T(t - 1) + ex->get_T(t - 1))
								/ (h * h * ex->get_mu() * (ex->get_mu() + 1)));
						else
							sum += h_t * ((nodes[i][j - 1]->get_T(t - 1) - 2 * nodes[i][j]->get_T(t - 1) + nodes[i][j + 1]->get_T(t - 1))
								/ (h * h));
						if (ey && ey->y < (nodes[i][j]->y) * h)
							sum += h_t * (2 * ((ey->get_mu()) * (nodes[i + 1][j]->get_T(t - 1)) - ((ey->get_mu()) + 1) * (nodes[i][j]->get_T(t - 1)) + (ey->get_T(t - 1)))
								/ (h * h * (ey->get_mu()) * ((ey->get_mu()) + 1)));
						else if (ey)
							sum += h_t * (2 * ((ey->get_mu()) * (nodes[i - 1][j]->get_T(t - 1)) - ((ey->get_mu()) + 1) * (nodes[i][j]->get_T(t - 1)) + (ey->get_T(t - 1)))
								/ (h * h * (ey->get_mu()) * ((ey->get_mu()) + 1)));
						else
							sum += h_t * ((nodes[i - 1][j]->get_T(t - 1) - 2 * nodes[i][j]->get_T(t - 1) + nodes[i + 1][j]->get_T(t - 1))
								/ (h * h));
						nodes[i][j]->set_T(t, sum);
					}
			for (int i = 0; i < edges.size(); i++)
				edges[i]->calculate(t);
		}
	}

	void calculate_implicit()  // Вычисление значений температур в узлах не явным методом с помощью метода расщепления
	{
		double mu;
		double* difw = new double(0);
		double* difr = new double(0);
		double** ws;
		ws = new double* [m];
		int size_matr;
		for (int i = 0; i < edges.size(); i++)
			edges[i]->calculate(0);
		for (int t = 1; t <= sim_t; t++)
		{
			for (int i = 0; i < m; i++)
			{
				ws[i] = new double[n];
				size_matr = 0;
				for (int j = 0; j < n; j++)
					if (nodes[i][j] && nodes[i][j]->is_active())
						size_matr++;
				Eigen::SparseMatrix<double> sm(size_matr, size_matr);
				Eigen::VectorXd res(size_matr), x(size_matr);
				size_matr = 0;
				for (int j = 0; j < n; j++)
				{
					if (nodes[i][j] && nodes[i][j]->is_active())
					{
						mu = 1;
						*difw = 0;
						*difr = 0;
						EdgeNode* ex = nodes[i][j]->get_edge_x();
						if (ex)
						{
							ex->get_edge_dif(difw, difr);
							mu = ex->get_mu();
						}
						if (!ex || (ex && ex->x < nodes[i][j]->x))
						{
							sm.insert(size_matr, size_matr + 1) = 2 * h_t / (h * h * (1 + mu));
						}
						if (!ex || (ex && ex->x >= nodes[i][j]->x))
						{
							sm.insert(size_matr, size_matr - 1) = 2 * h_t / (h * h * (1 + mu));
						}
						sm.insert(size_matr, size_matr) = -(2 * h_t / (mu * h * h) + 1) + *difw;
						res(size_matr) = -nodes[i][j]->get_T(t - 1) + *difr;
						size_matr++;
					}
				}
				if (size_matr > 0)
				{
					Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
					solver.compute(sm);
					x = solver.solve(res);
					size_matr = 0;
					for (int j = 0; j < n; j++)
					{
						if (nodes[i][j] && nodes[i][j]->is_active())
						{
							ws[i][j] = x(size_matr);
							size_matr++;
						}
					}
				}
			}

			for (int j = 0; j < n; j++)
			{
				size_matr = 0;
				for (int i = 0; i < m; i++)
					if (nodes[i][j] && nodes[i][j]->is_active())
						size_matr++;
				Eigen::SparseMatrix<double> sm(size_matr, size_matr);
				Eigen::VectorXd res(size_matr), x(size_matr);
				size_matr = 0;
				for (int i = 0; i < m; i++)
				{
					if (nodes[i][j] && nodes[i][j]->is_active())
					{
						mu = 1;
						*difw = 0;
						*difr = 0;
						EdgeNode* ey = nodes[i][j]->get_edge_y();
						if (ey)
						{
							ey->get_edge_dif(difw, difr);
							mu = ey->get_mu();
						}
						if (!ey || (ey && ey->y < nodes[i][j]->y))
						{
							sm.insert(size_matr, size_matr + 1) = 2 * h_t / (h * h * (1 + mu));
						}
						if (!ey || (ey && ey->y >= nodes[i][j]->y))
						{
							sm.insert(size_matr, size_matr - 1) = 2 * h_t / (h * h * (1 + mu));
						}
						sm.insert(size_matr, size_matr) = -(2 * h_t / (mu * h * h) + 1) + *difw;
						res(size_matr) = -ws[i][j] + *difr;
						size_matr++;
					}
				}
				if (size_matr > 0)
				{
					Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
					solver.compute(sm);
					x = solver.solve(res);
					size_matr = 0;
					for (int i = 0; i < m; i++)
					{
						if (nodes[i][j] && nodes[i][j]->is_active())
						{
							nodes[i][j]->set_T(t, x(size_matr));
							size_matr++;
						}
					}
				}
			}
			for (int i = 0; i < edges.size(); i++)
				edges[i]->calculate(t);
		}
	}

	void T_out(string name)  // Вывод результатов в файл
	{
		ofstream fout(name);
		
		for (int t = 0; t <= sim_t; t++)
			for (int i = 0; i < m; i++)
				for (int j = 0; j < n; j++)
					if (nodes[i][j] && nodes[i][j]->is_active())
						fout << t << " " << j * h << " " << i * h << " " << nodes[i][j]->get_T(t) << '\n';
					
		fout.close();
	}

	void show_bin()
	{
		cout << "show_bin()\n";
		for (int i = m - 1; i > -1; i--)
		{
			for (int j = 0; j < n; j++)
				cout << bin_mask[i][j] << ' ';
			cout << '\n';
		}
	}
	void show_form()
	{
		cout << "show_form()\n";
		for (int i = m - 1; i > -1; i--)
		{
			for (int j = 0; j < n; j++)
				if (nodes[i][j] && (nodes[i][j])->get_parent())
				{
					cout << (nodes[i][j]->get_parent())->get_form() << ' ';
				}
				else
					cout << "- ";
			cout << '\n';
		}
	}
	void show_active()
	{
		cout << "show_active()\n";
		for (int i = m - 1; i > -1; i--)
		{
			if (i < 10)
				cout << i << "  ";
			else
				cout << i << " ";
			for (int j = 0; j < n; j++)
				if (nodes[i][j])
				{
					cout << nodes[i][j]->is_active() << ' ';
				}
				else
					cout << "- ";
			cout << '\n';
		}
		cout << "   ";
		for (int j = 0; j < n; j++)
			if (j < 10)
				cout << j << " ";
			else
				cout << j;
		cout << '\n';
	}
	void show_edges()
	{
		cout << "show_edges()\n";
		cout << edges.size() << "!" << endl;
		for (int i = 0; i < edges.size(); i++)
			cout << edges[i]->get_mu() << " ";
		cout << endl;
	}
	
	void show_edge_cond()
	{
		cout << "show_edge_cond()\n";
		for (int i = m - 1; i > -1; i--)
		{
			for (int j = 0; j < n; j++)
				if (nodes[i][j] && nodes[i][j]->get_edge_x() && nodes[i][j]->get_edge_y())
					cout << ((nodes[i][j]->get_edge_x())->get_cond()).get_type() << ((nodes[i][j]->get_edge_y())->get_cond()).get_type() << " ";
				else if (nodes[i][j] && nodes[i][j]->get_edge_x())
					cout << ((nodes[i][j]->get_edge_x())->get_cond()).get_type() << "  ";
				else if (nodes[i][j] && nodes[i][j]->get_edge_y())
					cout << ((nodes[i][j]->get_edge_y())->get_cond()).get_type() << "  ";
				else if (nodes[i][j])
					cout << "00 ";
				else
					cout << "-- ";
			cout << endl;
		}
	}
	void show_edge_T(int t)
	{
		cout << "show_edge_T(" << t << ")\n";
		for (int i = m - 1; i > -1; i--)
		{
			for (int j = 0; j < n; j++)
				if (nodes[i][j] && nodes[i][j]->get_edge_x() && nodes[i][j]->get_edge_y())
					cout << ((nodes[i][j]->get_edge_x())->get_T(t)) << ((nodes[i][j]->get_edge_y())->get_T(t)) << " ";
				else if (nodes[i][j] && nodes[i][j]->get_edge_x())
					cout << ((nodes[i][j]->get_edge_x())->get_T(t)) << "  ";
				else if (nodes[i][j] && nodes[i][j]->get_edge_y())
					cout << ((nodes[i][j]->get_edge_y())->get_T(t)) << "  ";
				else if (nodes[i][j])
					cout << "00 ";
				else
					cout << "-- ";
			cout << endl;
		}
	}
	
	void show_T(int t)
	{
		cout << "show_T(" << t << ")\n";
		for (int i = m - 1; i > -1; i--)
		{
			for (int j = 0; j < n; j++)
				if (nodes[i][j] && nodes[i][j]->is_active())
					printf("%3.0f ", nodes[i][j]->get_T(t));
				else
					printf("--- ");
			cout << endl;
		}
	}
	void show_edge_mu()
	{
		cout << "show_edge_mu()\n";
		for (int i = m - 1; i > -1; i--)
		{
			for (int j = 0; j < n; j++)
				if (nodes[i][j] && nodes[i][j]->get_edge_x() && nodes[i][j]->get_edge_y())
					cout << ((nodes[i][j]->get_edge_x())->get_mu()) << ((nodes[i][j]->get_edge_y())->get_mu()) << " ";
				else if (nodes[i][j] && nodes[i][j]->get_edge_x())
					cout << ((nodes[i][j]->get_edge_x())->get_mu()) << "  ";
				else if (nodes[i][j] && nodes[i][j]->get_edge_y())
					cout << ((nodes[i][j]->get_edge_y())->get_mu()) << "  ";
				else if (nodes[i][j])
					cout << "00 ";
				else
					cout << "-- ";
			cout << endl;
		}
	}
};


int main()
{
	int h = 10;

	cout << "Input h:\n";
	cin >> h;
	if (h < 0)
		return 0;

	Object obj(500, 400, h, 100, 1);

	Rectangle r1(0, 0, 350, 400, EdgeCondition(3, 100)), r2(340, 0, 500, 250, EdgeCondition(3, 100));
	Circle c1(350, 250, 150, EdgeCondition(3, 200));


	Point pin(355, 255);  // Центр отверстия

	/// EdgeCondition определяет, какое из условий применяется, нумерация как в задании
	/// EdgeCondition(3, 100) => T = 100; EdgeCondition(4) => grad(T)=0; EdgeCondition(5) => grad(T)=T

	//Rectangle rin(pin.x - 50, pin.y - 50, pin.x + 50, pin.y + 50, EdgeCondition(5));  /// Отверстие - круг
	Circle ci(pin.x, pin.y, 50, EdgeCondition(4));  /// Отверстие - квадрат

	/// Граничные условия, r, l, t, b определяют сторону
	Rectangle rr(500, 0, 500, 250, EdgeCondition(5)), rl(0, 0, 0, 400, EdgeCondition(3, 100));
	Rectangle rt(0, 400, 350, 400, EdgeCondition(5)), rb(0, 0, 500, 0, EdgeCondition(3, 100));
	
	obj.add_form(&r1, 0, 1);
	obj.add_form(&r2, 0, 1);
	obj.add_form(&c1, 0, 1);
	
	obj.add_form(&ci, 0, 0);
	obj.add_form(&rr, 0, 0);
	obj.add_form(&rl, 0, 0);
	obj.add_form(&rt, 0, 0);
	obj.add_form(&rb, 0, 0);

	/*obj.show_bin();
	obj.show_active();*/
	//obj.show_form();

	obj.compile_edges();

	/*obj.show_bin();
	obj.show_active();
	
	obj.show_edges();
	obj.show_edge_cond();
	obj.show_edge_mu();*/

	string file_ex = "explicit", file_im = "implicit", h_str;
	int h_copy = h;
	while (h_copy > 0)
	{
		h_str += '0' + h_copy % 10;
		h_copy /= 10;
	}
	h_str = string(h_str.rbegin(), h_str.rend());
	file_ex += h_str;
	file_im += h_str;
	file_ex += ".txt";
	file_im += ".txt";

	cout << "Explicit method results:" << endl;
	obj.calculate_explicit();
	//obj.show_edge_T(100);
	obj.show_T(100);
	obj.T_out(file_ex);

	cout << "Implicit method results:" << endl;
	obj.calculate_implicit();
	//obj.show_edge_T(100);
	obj.show_T(100);
	obj.T_out(file_im);
}