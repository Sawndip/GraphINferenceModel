string int_to_str(int i){
   string s;
   stringstream out;
   out << i;
   s = out.str();
   return s;
}


string bool_to_str (bool b)
{
   return b ? "true" : "false";
}


string double_to_str(double value) {
   stringstream sstr;
   sstr << value;
   return sstr.str();
}


float string_to_float(string str) {
//   ifastream<basic_formatters, string_reader> myString(&str);
//   float value;
//   myString >> value;
//   return value;
  return atof(str.c_str());
}

double string_to_double(string str) {
    //   ifastream<basic_formatters, string_reader> myString(&str);
    //   float value;
    //   myString >> value;
    //   return value;
    return atof(str.c_str());
}
