"""A minimal stand-in for the browser's `pyscript` module, used only by tests.

Like real PyScript, @when wraps a handler so it must be called with an event.
"""
import inspect


class Element:
    def __init__(self, id_, value=""):
        self.id, self.value = id_, value
        self.innerText = self.innerHTML = self.className = self.src = ""
        self.disabled = self.readOnly = False
        self.style = type("Style", (), {"display": ""})()


DEFAULTS = dict(L="1.0", h="0.1", t="0.01", p="100", load="bending", preset="steel_a36",
                element="Q4", nx="40", ny="4", field="von_mises")


class Document:
    def __init__(self):
        self.elements = {}

    def getElementById(self, id_):
        if id_ not in self.elements:
            self.elements[id_] = Element(id_, DEFAULTS.get(id_, ""))
        return self.elements[id_]


document = Document()
handlers = {}


def when(event_type, selector):
    def decorate(func):
        if inspect.iscoroutinefunction(func):
            async def wrapper(event):
                return await func(event)
        else:
            def wrapper(event):
                return func(event)
        handlers[selector] = wrapper
        return wrapper
    return decorate
