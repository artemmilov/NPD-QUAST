def parse_from_mgf(s):
    res = []
    cur = ''
    in_brackets = False
    for a in s:
        if (a == ',') and (not in_brackets):
            res.append(cur)
            cur = ''
        elif a == '\"':
            in_brackets = not in_brackets
        else:
            cur += a
    if cur != '':
        res.append(cur)
    return res


def get_name(conn, user_id):
    c = conn.cursor()
    c.execute("SELECT name FROM userList WHERE userID = ?", (user_id,))
    result = c.fetchone()
    if result:
        return result[0]
