// THE BINDINGS OF TRIADTITANS/BLA, READY TO BE #INCLUDED INTO A PYBIND11_MODULE MACRO

m.doc() = "Basic linear algebra module"; // optional module docstring

py::class_<Vector<double>> (m, "Vector", py::module_local())
    .def(py::init<size_t>(),
        py::arg("size"), "create vector of given size")
    .def("__len__", &Vector<double>::Size,
        "return size of vector")
    
    .def("__setitem__", [](Vector<double> & self, int i, double v) {
    if (i < 0) i += self.Size();
    if (i < 0 || i >= self.Size()) throw py::index_error("vector index out of range");
    self(i) = v;
    })
    .def("__getitem__", [](Vector<double> & self, int i) { return self(i); })
    
    .def("__setitem__", [](Vector<double> & self, py::slice inds, double val)
    {
    size_t start, stop, step, n;
    if (!inds.compute(self.Size(), &start, &stop, &step, &n))
        throw py::error_already_set();
    self.Range(start, stop+1).Slice(0,step) = val;
    })
    
    .def("__add__", [](Vector<double> & self, Vector<double> & other)
    { return Vector<double> (self+other); })

    .def("__rmul__", [](Vector<double> & self, double scal)
    { return Vector<double> (scal*self); })
    
    .def("__str__", [](const Vector<double> & self)
    {
    std::stringstream str;
    str << self;
    return str.str();
    })

    .def(py::pickle(
    [](Vector<double> & self) { // __getstate__
        /* return a tuple that fully encodes the state of the object */
        return py::make_tuple(self.Size(),
                            py::bytes((char*)(void*)&self(0), self.Size()*sizeof(double)));
    },
    [](py::tuple t) { // __setstate__
        if (t.size() != 2)
        throw std::runtime_error("should be a 2-tuple!");

        Vector<double> v(t[0].cast<size_t>());
        py::bytes mem = t[1].cast<py::bytes>();
        std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()), v.Size()*sizeof(double));
        return v;
    }))
;

py::class_<Vec<3>> (m, "Vec3", py::module_local())
    .def(py::init<>(), "create vector of fixed size 3")
    .def("__len__", &Vec<3>::Size,
        "return size of vector")
    
    .def("__setitem__", [](Vec<3> & self, int i, double v) {
    if (i < 0) i += self.Size();
    if (i < 0 || i >= self.Size()) throw py::index_error("vector index out of range");
    self(i) = v;
    })
    .def("__getitem__", [](Vec<3> & self, int i) { return self(i); })
    
    .def("__setitem__", [](Vec<3> & self, int i, double val)
    {
    self(i) = val;
    })
    
    .def("__str__", [](const Vec<3> & self)
    {
    std::stringstream str;
    str << self;
    return str.str();
    })
;

py::class_<Matrix<double,Ordering::RowMajor>> (m, "Matrix", py::buffer_protocol(), py::module_local())
    .def(py::init<size_t,size_t>(),
        py::arg("height"), py::arg("width"), "create matrix of given size")
    .def("__setitem__", [](Matrix<double,Ordering::RowMajor> & self, std::tuple<int, int> shape, double v) {
    auto [row,col] = shape;
    if (row < 0) row += self.Height();
    if (col < 0) col += self.Width();
    // row, col now mutated

    if (row < 0 || row >= self.Height()) throw py::index_error("matrix row out of range");
    if (col < 0 || col >= self.Width()) throw py::index_error("matrix col out of range");
    self(row, col) = v;
    })
    .def("__getitem__", [](Matrix<double,Ordering::RowMajor> & self, std::tuple<int,int> shape) {
    auto [row, col] = shape;
    return self(row, col);
    })
    .def("inverse", [](Matrix<double,Ordering::RowMajor> & self) {
    return inverse(self);
    })

    .def_buffer([](Matrix<double,Ordering::RowMajor> &m) -> py::buffer_info {
    std::cout << "Using buffer protocol ...";
    return py::buffer_info(
        &m(0,0),                               /* Pointer to buffer */
        sizeof(double),                          /* Size of one scalar */
        py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
        2,                                      /* Number of dimensions */
        { m.Height(), m.Width() },                 /* Buffer dimensions */
        { sizeof(double) * m.Dist(),             /* Strides (in bytes) for each index */
            sizeof(double) }
    );
    })
    
    // We don't know yet what slicing should do with matrices.
    // .def("__setitem__", [](Vector<double> & self, py::slice inds, double val)
    // {
    //   size_t start, stop, step, n;
    //   if (!inds.compute(self.Size(), &start, &stop, &step, &n))
    //     throw py::error_already_set();
    //   self.Range(start, stop).Slice(0,step) = val;
    // })
    
    .def("__add__", [](Matrix<double> & self, Matrix<double> & other)
    { return Matrix<double> (self+other); })

    .def_property_readonly("shape",
    [](const Matrix<double, Ordering::RowMajor>& self) {
        return std::tuple(self.Height(), self.Width());
    })

    .def("__rmul__", [](Matrix<double> & self, double scal)
    { return Matrix<double> (scal*self); })
    
    .def("__mul__", [](Matrix<double> & self, Matrix<double> other)
    { return Matrix<double> (self*other); })

    .def("__mul__", [](Matrix<double> & self, Vector<double> other)
    { return Vector<double> (self*other); })

    .def("__str__", [](const Matrix<double,Ordering::RowMajor> & self)
    {
    std::stringstream str;
    str << self;
    return str.str();
    })

    .def(py::pickle(
    [](Matrix<double> & self) { // __getstate__
        /* return a tuple that fully encodes the state of the object */
        return py::make_tuple(self.Height(),self.Width(),
                            py::bytes((char*)(void*)&self(0,0), self.Height()*self.Width()*sizeof(double)));
    },
    [](py::tuple t) { // __setstate__
        if (t.size() != 3)
        throw std::runtime_error("should be a 3-tuple!");

        Matrix<double> m(t[0].cast<size_t>(),t[1].cast<size_t>());
        py::bytes mem = t[2].cast<py::bytes>();
        std::memcpy(&m(0,0), PYBIND11_BYTES_AS_STRING(mem.ptr()), m.Height()*m.Width()*sizeof(double));
        return m;
    }))
;