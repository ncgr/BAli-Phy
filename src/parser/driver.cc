#include "driver.hh"
#include "parser.hh"

void driver::pop_context()
{
    contexts.pop_back();
}

boost::optional<LayoutContext> driver::get_context()
{
    if (contexts.empty())
	return boost::none;
    return contexts.back();
}

void driver::push_context(const boost::optional<LayoutContext>& lc)
{
    contexts.push_back(lc);
}

void driver::push_context(const LayoutContext& lc)
{
    contexts.push_back(lc);
}

void driver::push_context()
{
    contexts.push_back({});
}

LayoutContext driver::get_offside(const yy::parser::location_type& loc)
{
    return {0,true};
}

driver::driver ()
    : trace_parsing (false), trace_scanning (false)
{
}

int
driver::parse (const std::string &f)
{
  file = f;
  location.initialize (&file);
  scan_begin ();
  yy::parser parser (*this);
  parser.set_debug_level (trace_parsing);
  int res = parser.parse ();
  scan_end ();
  return res;
}

